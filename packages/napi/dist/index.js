"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.verify_with_file_path = exports.prove_with_file_path = exports.verify = exports.prove = void 0;
var _a = require("../index.node"), _prove = _a.prove, _verify = _a.verify, _prove_with_file_path = _a.prove_with_file_path, _verify_with_file_path = _a.verify_with_file_path;
function prove(r1cs_file_path, witness_json_path, proof_json_path) {
    return _prove(r1cs_file_path, witness_json_path, proof_json_path);
}
exports.prove = prove;
function verify(r1cs_file_path, witness_json_path, proof_json_path) {
    return _verify(r1cs_file_path, witness_json_path, proof_json_path);
}
exports.verify = verify;
function prove_with_file_path(r1cs_file_path, witness_json_path, proof_json_path) {
    return _prove_with_file_path(r1cs_file_path, witness_json_path, proof_json_path);
}
exports.prove_with_file_path = prove_with_file_path;
function verify_with_file_path(r1cs_file_path, witness_json_path, proof_json_path) {
    return _verify_with_file_path(r1cs_file_path, witness_json_path, proof_json_path);
}
exports.verify_with_file_path = verify_with_file_path;
//# sourceMappingURL=index.js.map