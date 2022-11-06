<< xDecomp`

$DefInfoQ = False;

DefConstantSymbol[dimx];
DefParameter /@ {r, t};
DefScalarFunction /@ {h, f};
DefManifold[Mf, dimx + 2, IndexRangeNS[Mf`A, Mf`H]];
DefManifold[Mx, dimx, IndexRangeNS[Mx`a, Mx`h]];
DefMetric[NoSignDet, metricMf[-Mf`A, -Mf`B], CDMf];
DefMetric[NoSignDet, metricMx[-Mx`a, -Mx`b], CDMx];
DefGBasis[decomp, {{et[Mf`A], edt[-Mf`A], t}, {er[Mf`A], edr[-Mf`A], r}}, {eX[Mx`a, -Mf`A]}];

MUnit`BeginTestSection["Decomposing SSS metric ansatze"];

SetHeldMetric[covd1, metricMf,
    CreateGCTensor[{
        {-t, -t} -> -h[r],
        {-r, -r} -> 1 / f[r],
        {-Mx`a, -Mx`b} -> r^2 metricMx[-Mx`a, -Mx`b]
    }, -{decomp, decomp}]
,
    CreateGCTensor[{
        {t, t} -> -1 / h[r],
        {r, r} -> f[r],
        {Mx`a, Mx`b} -> 1 / r^2 metricMx[Mx`a, Mx`b]
    }, {decomp, decomp}]
];
AddCurvatureTensorsToHolder[covd1, decomp, ChristoffelCDMf];
VerificationTest[
    CachedGCTensor[covd1, RicciScalarCDMf, {}][],
    (1/(2 r^2 h[r]^2))(h[r] (2 h[
       r] (RicciScalarCDMx[] - dimx r f'[r]) -
     r^2 f'[r] h'[r]) +
  f[r] (-2 (-1 + dimx) dimx h[r]^2 + r^2 h'[r]^2 -
     2 r h[r] (dimx h'[r] + r h''[r])))
];

VerificationTest[
    (RiemannCDMf[-Mf`A, Mf`C, -Mf`B, Mf`D]
    * RiemannCDMf[-Mf`C, Mf`E, -Mf`D, Mf`F]
    * RiemannCDMf[-Mf`E, Mf`A, -Mf`F, Mf`B]
    /. GetAllHeldTensorRules[covd1]
    // ContractGCTensors[covd1]
    // NoScalar
    // ContractMetric
    // Simplification)
    - (1/(8 r^6 h[r]^3) (2 h[r]^3 (4 RiemannCDMx[-Mx`a, Mx`e, -Mx`c, Mx`f]RiemannCDMx[Mx`a, Mx`b, Mx`c, Mx`d]
       RiemannCDMx[-Mx`b, -Mx`e, -Mx`d, -Mx`f] +
      3 r^2 RicciScalarCDMx[] f'[r]^2) -
   2 (-1 + dimx) dimx f[r]^3 h[r] (4 (-2 + dimx) h[r]^2 + 3 r^2 h'[r]^2) - 
   3 f[r] (2 h[r]^3 (4  RicciCDMx[-Mx`a, -Mx`b]RicciCDMx[Mx`a, Mx`b] - 
         4  RiemannCDMx[-Mx`a, -Mx`c, -Mx`b, -Mx`d] RiemannCDMx[Mx`a, Mx`b, Mx`c, Mx`d] + (-1 + dimx) dimx r^2 f'[r]^2) + 
      dimx r^4 h[r] f'[r]^2 h'[r]^2) + 
   3 f[r]^2 (8 (-2 + dimx) h[r]^3 RicciScalarCDMx[] + 
      dimx r^4 f'[r] h'[r]^3 + 
      2 r^2 h[r] h'[
        r] (RicciScalarCDMx[] h'[r] - 
         dimx r^2 f'[r] h''[r])))
    ) // Simplification
,
    0
];

MUnit`EndTestSection[];

UndefTensor /@ {et, edt, er, edr, eX};
Undef /@ VisitorsOf@metricMf;
Undef /@ VisitorsOf@metricMx;
UndefMetric /@ {metricMx, metricMf};
Undef /@ VisitorsOf@Mx;
Undef /@ VisitorsOf@Mf;
UndefManifold /@ {Mx, Mf};
UndefConsantSymbol[dimx];
UndefParameter /@ {r, t};