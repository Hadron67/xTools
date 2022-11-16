Needs["xTools`Misc`", "Misc.wl"];

MUnit`BeginTestSection["FVariation"];

VerificationTest[FVariation[2 f[x]g'[x]], 2 FVariation[f[x]] g'[x] + 2 f[x] FVariation[g'[x]]];
VerificationTest[FVariation[2 f[x]g'[x] g[x], 1, ConstantFunctions -> {g}], 2 FVariation[f[x], 1, ConstantFunctions -> {g}] g'[x] g[x]];

MUnit`EndTestSection[];