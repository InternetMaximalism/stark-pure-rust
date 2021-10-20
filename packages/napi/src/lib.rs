use neon::prelude::*;
use r1cs_stark::run::{
    prove_with_file_path as _prove_with_file_path, verify_with_file_path as _verify_with_file_path,
};

fn prove_with_file_path(mut cx: FunctionContext) -> JsResult<JsUndefined> {
    let r1cs_file_path = cx.argument::<JsString>(0)?.value(&mut cx);
    let witness_file_path = cx.argument::<JsString>(1)?.value(&mut cx);
    let proof_json_path = cx.argument::<JsString>(2)?.value(&mut cx);
    _prove_with_file_path(&r1cs_file_path, &witness_file_path, &proof_json_path).unwrap();
    Ok(cx.undefined())
}

// "../r1cs-stark/tests/mul_bn128.r1cs", "../r1cs-stark/tests/mul_bn128_wtns_valid.json", "../r1cs-stark/tests/mul_bn128_proof.json"
fn verify_with_file_path(mut cx: FunctionContext) -> JsResult<JsBoolean> {
    let r1cs_file_path = cx.argument::<JsString>(0)?.value(&mut cx);
    let witness_file_path = cx.argument::<JsString>(1)?.value(&mut cx);
    let proof_json_path = cx.argument::<JsString>(2)?.value(&mut cx);
    _verify_with_file_path(&r1cs_file_path, &witness_file_path, &proof_json_path).unwrap();
    Ok(cx.boolean(true))
}

#[neon::main]
fn main(mut cx: ModuleContext) -> NeonResult<()> {
    cx.export_function("prove", prove_with_file_path)?;
    cx.export_function("verify", verify_with_file_path)?;
    cx.export_function("prove_with_file_path", prove_with_file_path)?;
    cx.export_function("verify_with_file_path", verify_with_file_path)?;
    Ok(())
}
