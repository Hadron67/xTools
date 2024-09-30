Needs["xTools`Misc`", "Misc.wl"];

VerificationTest[FVariation[2 f[x]g'[x]], 2 FVariation[f[x]] g'[x] + 2 f[x] FVariation[g'[x]]];
VerificationTest[FVariation[2 f[x]g'[x] g[x], 1, ConstantFunctions -> {g}], 2 FVariation[f[x], 1, ConstantFunctions -> {g}] g'[x] g[x]];

With[{
    case = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}
},
    VerificationTest[SparseRowReduce@case, RowReduce@case];
];

TestSparseRowReduce[mat_] := VerificationTest[RowReduce@mat, SparseRowReduce@mat];

Table[TestSparseRowReduce@Table[RandomInteger[100], 10, 10], 10];
