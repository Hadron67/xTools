<< xDecomp`

$DefInfoQ = False;

MUnit`BeginTestSection["GCTensor"];

DefConstantSymbol[dimx];
DefManifold[Mf, dimx + 2, {Mf`A, Mf`B, Mf`C, Mf`D, Mf`E, Mf`F, Mf`G, Mf`H}];
DefManifold[Mx, dimx, {Mx`a, Mx`b, Mx`c, Mx`d, Mx`e, Mx`f, Mx`g, Mx`h}];
DefMetric[NoSignDet, metricMf[-Mf`A, -Mf`B], CDMf];


MUnit`EndTestSection[];