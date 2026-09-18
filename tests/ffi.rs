//! In-process smoke tests for the C ABI: flat arrays in, structs out.

use max_inscribed_circle::ffi::{
    mic_options_default, mic_solve, mic_solve_multipolygon, MicOptionsC, MicPolyInput, MicResultC,
};
use std::mem::MaybeUninit;

fn square() -> Vec<f64> {
    vec![
        0.0, 0.0, //
        10.0, 0.0, //
        10.0, 10.0, //
        0.0, 10.0, //
        0.0, 0.0,
    ]
}

#[test]
fn ffi_mic_square() {
    let sq = square();
    unsafe {
        let mut out = MaybeUninit::<MicResultC>::uninit();
        let rc = mic_solve(
            sq.as_ptr(),
            10,
            std::ptr::null(),
            0,
            std::ptr::null(),
            out.as_mut_ptr(),
        );
        assert_eq!(rc, 0);
        let r = out.assume_init();
        assert!((r.radius - 5.0).abs() < 1e-6, "radius={}", r.radius);
        assert!((r.center_x - 5.0).abs() < 1e-6);
        assert!((r.center_y - 5.0).abs() < 1e-6);
        // A square is answered by the analytical fast path, so engine code 0.
        assert_eq!(r.used_engine, 0);
    }
}

#[test]
fn ffi_mic_multipolygon_keeps_best() {
    let big = square();
    let small = vec![
        20.0, 20.0, //
        22.0, 20.0, //
        22.0, 22.0, //
        20.0, 22.0, //
        20.0, 20.0,
    ];
    let inputs = [
        MicPolyInput {
            coords: small.as_ptr(),
            coords_len: 10,
            holes: std::ptr::null(),
            holes_len: 0,
        },
        MicPolyInput {
            coords: big.as_ptr(),
            coords_len: 10,
            holes: std::ptr::null(),
            holes_len: 0,
        },
    ];
    unsafe {
        let mut out = MaybeUninit::<MicResultC>::uninit();
        let rc = mic_solve_multipolygon(inputs.as_ptr(), 2, std::ptr::null(), out.as_mut_ptr());
        assert_eq!(rc, 0);
        let r = out.assume_init();
        assert!((r.radius - 5.0).abs() < 1e-6, "radius={}", r.radius);
        assert_eq!(r.component_index, 1);
    }
}

#[test]
fn ffi_mic_rejects_garbage() {
    unsafe {
        let mut out = MaybeUninit::<MicResultC>::uninit();
        let rc = mic_solve(
            std::ptr::null(),
            0,
            std::ptr::null(),
            0,
            std::ptr::null(),
            out.as_mut_ptr(),
        );
        assert_eq!(rc, -1);
    }
}

#[test]
fn ffi_mic_honours_explicit_options() {
    // Explicit options: FastF64 robustness, solved through the native pipeline.
    let opts = MicOptionsC { robust_mode: 0 };
    let sq = square();
    unsafe {
        let mut out = MaybeUninit::<MicResultC>::uninit();
        let rc = mic_solve(
            sq.as_ptr(),
            10,
            std::ptr::null(),
            0,
            &opts,
            out.as_mut_ptr(),
        );
        assert_eq!(rc, 0);
        let r = out.assume_init();
        assert!((r.radius - 5.0).abs() < 1e-6, "radius={}", r.radius);
        assert_eq!(r.used_engine, 0);
    }
    assert_eq!(mic_options_default().robust_mode, 1);
}
