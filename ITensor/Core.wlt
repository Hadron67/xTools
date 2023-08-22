<< xAct`xPerm`
<< xTools`ITensor`Core`

SlotsOfIObject[t0, inds___] ^:= ConstantArray[dimf, Length@{inds}];
SlotsOfIObject[Riemann, inds___] ^:= With[{len = Length@{inds}}, Join[ConstantArray[ITensorToScalarMapSlot, len - 4], ConstantArray[dimf, 4]]];
SymmetryGroupOfIObject[Riemann, _, _, _, __] ^:= StrongGenSet[Range@4, GenSet[-xAct`xPerm`Cycles@{1, 2}, -xAct`xPerm`Cycles@{3, 4}, xAct`xPerm`Cycles[{1, 3}, {2, 4}]]]


MUnit`BeginTestSection["FindIndices"];

VerificationTest[TensorMapGroups[{ITensorToTensorMapSlot[1], ITensorToTensorMapSlot[1], ITensorToTensorMapSlot[2], ITensorToTensorMapSlot@None, ITensorToTensorMapSlot[2], ITensorToTensorMapSlot@None}], {{4}, {6}, {1, 2}, {3, 5}}];
VerificationTest[FindFreeIndices@t0[a, b, c, -c], {a, b}];
VerificationTest[FindFreeIndices[t0[a, b, c, -c] + t0[a, b]], {a, b}];
VerificationTest[FindFreeAndDummyIndices[t0[a, b, c, -c, d]t0[-d] + 3 t0[a, b, c]t0[-c]], {{a, b}, {c, d}}];
VerificationTest[FindFreeIndices[t0[c] ^ 2], {c}];
VerificationTest[FindFreeIndices[2 Riemann[e, -e, f, -f] Riemann[2, a, b, c, d]Riemann[1, -a, -b, -c , -d] t0[h]], {h}];

MUnit`EndTestSection[];

MUnit`BeginTestSection["FindIndices"];

VerificationTest[
    ISort[2 Riemann[e, -e, f, -f] Riemann[2, a, b, c, d]Riemann[1, -a, -b, -c , -d] t0[h]]
,
    ISortedIObject[Times, {
        ISortedIObject[Riemann, {e, -e, f, -f}, 8, 2, 0],
        ISortedIObject[Riemann, {1, -a, -b, -c, -d}, 8, 0, 4],
        ISortedIObject[Riemann, {2, a, b, c, d}, 8, 0, 4],
        ISortedIObject[t0, {h}, 1, 0, 1],
        ISortedOther[2]
    }, 1, 6, 1]
];

MUnit`EndTestSection[];
