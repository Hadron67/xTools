<< xDecomp`

$DefInfoQ = False;

MUnit`BeginTestSection["GCTensor"];

DefConstantSymbol[dimx];
DefManifold[Mf, dimx, {Mf`A, Mf`B, Mf`C, Mf`D, Mf`E, Mf`F, Mf`G, Mf`H}];
DefMetric[-1, metricMf[-Mf`A, -Mf`B], CDMf];

MUnit`EndTestSection[];