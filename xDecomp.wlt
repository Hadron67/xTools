Needs["xTools`xDecomp`", "xDecomp.wl"];

$DefInfoQ = False;

DefConstantSymbol[dimx];
DefParameter /@ {r, t};
r /: PD[_]@r = 0;
t /: PD[_]@t = 0;
DefScalarFunction /@ {h, f, Pn};
DefManifold[Mf, dimx + 2, IndexRangeNS[Mf`A, Mf`H]];
DefManifold[Mx, dimx, IndexRangeNS[Mx`a, Mx`h]];
DefMetric[NoSignDet, metricMf[-Mf`A, -Mf`B], CDMf];
DefMetric[NoSignDet, metricMx[-Mx`a, -Mx`b], CDMx];

DefGChart[decomp, TangentMf, {2, {TangentMx}}, PDGChartInfo[{t, r}, {TangentMx}]];
DefTensor[n0[Mf`A], {Mf, r}];
DefTensor[n1[Mf`A], {Mf, r}];
DefTensor[v0[Mx`a], {Mx, r}];
DefTensor[phi0[], {Mf, r, t}];

DefConstantSymbol[L0];
DefParameter[r2];
r2 /: PD[_]@r2 = 0;
DefManifold[Mr, dimx, IndexRangeNS[Mr`a, Mr`h]];
DefMetric[NoSignDet, metricMr[-Mr`a, -Mr`b], CDMr, OtherDependencies -> {r2}];
DefGChart[fgc, TangentMf, {1, {TangentMr}}, PDGChartInfo[{r2}, {TangentMr}]];

MUnit`BeginTestSection["Decomposing SSS metric ansatze"];

metricVal = GCTensorFromSparse[{
    {1, 1} -> -h[r],
    {2, 2} -> 1 / f[r],
    {3, 3} -> ETensor[r^2 metricMx[-Mx`a, -Mx`b]]
}, -{decomp, decomp}];

metricInvVal = GCTensorFromSparse[{
    {1, 1} -> -1 / h[r],
    {2, 2} -> f[r],
    {3, 3} -> ETensor[1 / r^2 metricMx[Mx`a, Mx`b]]
}, {decomp, decomp}];

detg = h[r] / f[r] r^(2dimx);

SetHeldMetric[holder1, metricMf, metricVal, metricInvVal];
AddCurvatureTensorsToHolder[holder1, decomp, ChristoffelCDMf];

VerificationTest[
    ContractMetric[CachedGCTensor[holder1, RicciScalarCDMf, {}][], {metricMx}] - (1/(2 r^2 h[r]^2))(h[r] (2 h[
       r] (RicciScalarCDMx[] - dimx r f'[r]) -
     r^2 f'[r] h'[r]) +
  f[r] (-2 (-1 + dimx) dimx h[r]^2 + r^2 h'[r]^2 -
     2 r h[r] (dimx h'[r] + r h''[r]))) // Simplify
,
    0
];

VerificationTest[
    (RiemannCDMf[-Mf`A, Mf`C, -Mf`B, Mf`D]
    * RiemannCDMf[-Mf`C, Mf`E, -Mf`D, Mf`F]
    * RiemannCDMf[-Mf`E, Mf`A, -Mf`F, Mf`B]
    // ContractGCTensors[holder1]
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

MUnit`BeginTestSection["Outer product"];

holder1 /: HeldGCTensor[holder1, n0] = GCTensorFromSparse[{{2} -> Sqrt[f[r]]}, {decomp}];
VerificationTest[
    ETensor[(n0[Mf`A] n0[Mf`B] // ContractGCTensors[holder1]) - GCTensorFromSparse[{{2, 2} -> f[r]}, {decomp, decomp}][Mf`A, Mf`B], {Mf`A, Mf`B}] // ScreenDollarIndices
,
    GCTensorFromSparse[{}, {decomp, decomp}]
];

VerificationTest[
    n0[Mf`A] n0[Mf`B] n0[-Mf`A] n0[-Mf`B] // ContractGCTensors[holder1]
,
    1
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["FG expansion"];

SetHeldMetric[fgcHolder, metricMf,
    GCTensorFromSparse[{
        {1, 1} -> L0^2/(4 r2^2),
        {2, 2} -> ETensor[1/r2 metricMr[-Mr`a, -Mr`b]]
    }, -{fgc, fgc}],
    GCTensorFromSparse[{
        {1, 1} -> (4 r2^2)/L0^2,
        {2, 2} -> ETensor[r2 metricMr[Mr`a, Mr`b]]
    }, {fgc, fgc}]
];
AddCurvatureTensorsToHolder[fgcHolder, fgc, ChristoffelCDMf];
GCTensorHolderDAUseMetricVB[fgcHolder] ^= None;

fgEEom = ETensor[
    RicciCDMf[-Mf`A, -Mf`B] + dimx / L0^2 metricMf[-Mf`A, -Mf`B]
    // ContractGCTensors[fgcHolder], {-Mf`A, -Mf`B}] // NoScalar // ToCanonical[#, UseMetricOnVBundle -> None] & // Simplify;
fgEEomExpected = GCTensor[
    {{1/4  metricMr[Mr`a, Mr`b] (metricMr[Mr`c, Mr`d]  ParamD[r2][metricMr[-Mr`a, -Mr`c]] ParamD[r2][metricMr[-Mr`b, -Mr`d]] - 2 ParamD[r2, r2][metricMr[-Mr`a, -Mr`b]]), 
    ETensor[1/2  metricMr[Mr`b, Mr`c]  (-CDMr[-Mr`a][
    ParamD[r2][metricMr[-Mr`b, -Mr`c]]] + CDMr[-Mr`c][
    ParamD[r2][metricMr[-Mr`a, -Mr`b]]]), {-Mr`a}]}, 
    {ETensor[1/2 metricMr[Mr`b, Mr`c]  (-CDMr[-Mr`a][
    ParamD[r2][
    metricMr[-Mr`b, -Mr`c]]] + CDMr[-Mr`c][
    ParamD[r2][
    metricMr[-Mr`a, -Mr`b]]]), {-Mr`a}], 
    ETensor[(1/(
        L0^2))(L0^2  RicciCDMr[-Mr`a, -Mr`b] + 
        2 r2  metricMr[Mr`c, Mr`d]  ParamD[r2][
    metricMr[-Mr`a, -Mr`c]] ParamD[r2][
    metricMr[-Mr`b, -Mr`d]] + 
        metricMr[-Mr`a, -Mr`b]   metricMr[Mr`c, Mr`d] 
            ParamD[r2][
    metricMr[-Mr`c, -Mr`d]] + ParamD[r2][
    metricMr[-Mr`a, -Mr`b]] (-2 + 
            dimx - r2  metricMr[Mr`c, Mr`d] 
            ParamD[r2][
    metricMr[-Mr`c, -Mr`d]]) - 
        2 r2 ParamD[r2, r2][
    metricMr[-Mr`a, -Mr`b]]), {-Mr`a, -Mr`b}]}}
, {-fgc, -fgc}];

VerificationTest[
    Simplify@ToCanonical[fgEEom - fgEEomExpected, UseMetricOnVBundle -> None]
,
    GCTensorFromSparse[{}, -{fgc, fgc}]
]

MUnit`EndTestSection[];

MUnit`BeginTestSection["ZeroGCTensorQ"];

VerificationTest[ZeroGCTensorQ@GCTensorFromSparse[{}, -{fgc, fgc, fgc}]];

MUnit`EndTestSection[];

MUnit`BeginTestSection["PostETensorContract"];

VerificationTest[
    GCTensorContractTwo[GCTensorFromSparse[{
        {3, 3} -> ETensor[metricMx[-Mx`a, -Mx`b] + metricMx[-Mx`b, -Mx`a]]
    }, -{decomp, decomp}], GCTensorFromSparse[{
        {3, 3} -> ETensor[metricMx[Mx`a, Mx`b] + metricMx[Mx`b, Mx`a]]
    }, {decomp, decomp}], {1}, {1}, PostETensorContract -> ContractMetric]
    - GCTensorFromSparse[{{3, 3} -> 4 ETensor[delta[-Mx`a, Mx`b]]}, {-decomp, decomp}] // ZeroGCTensorQ
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["Covariant derivatives"];

v1 = GCTensor[{
    h[t, r],
    f[t, r],
    ETensor[v0[Mx`a]]
}, {decomp}];

t1 = GCTensorFromSparse[{
    {1, 1} -> h[t, r],
    {2, 2} -> f[t, r],
    {2, 3} -> ETensor[f[t, r] v0[-Mx`a]],
    {3, 3} -> ETensor[v0[Mx`a] v0[-Mx`b]]
}, {decomp, -decomp}];

p0 = GCTensorFromSparse[{
    {1, 2, 1, 2} -> Pn[1, r],
    {1, 2, 2, 1} -> -Pn[1, r],
    {2, 1, 1, 2} -> -Pn[1, r],
    {2, 1, 2, 1} -> Pn[1, r],

    {1, 3, 1, 3} -> ETensor[Pn[2, r] metricMx[Mx`a, Mx`b]],
    {1, 3, 3, 1} -> ETensor[-Pn[2, r] metricMx[Mx`a, Mx`b]],
    {3, 1, 1, 3} -> ETensor[-Pn[2, r] metricMx[Mx`a, Mx`b]],
    {3, 1, 3, 1} -> ETensor[Pn[2, r] metricMx[Mx`a, Mx`b]],

    {2, 3, 2, 3} -> ETensor[Pn[3, r] metricMx[Mx`a, Mx`b]],
    {2, 3, 3, 2} -> ETensor[-Pn[3, r] metricMx[Mx`a, Mx`b]],
    {3, 2, 2, 3} -> ETensor[-Pn[3, r] metricMx[Mx`a, Mx`b]],
    {3, 2, 3, 2} -> ETensor[Pn[3, r] metricMx[Mx`a, Mx`b]],

    {3, 3, 3, 3} -> ETensor[Pn[4, r] (
        metricMx[Mx`a, Mx`c] metricMx[Mx`b, Mx`d]
        - metricMx[Mx`a, Mx`d] metricMx[Mx`b, Mx`c]
    )]
}, {decomp, decomp, decomp, decomp}];

covd = CovDGChart[-decomp, CachedGCTensor[holder1, ChristoffelCDMf, {1, -1, -1}]];
pd1 = PDGChart[-decomp];
chris = CachedGCTensor[holder1, ChristoffelCDMf, {1, -1, -1}];

VerificationTest[
    GCTensorPDDiv[v1, 1][] // NoScalar // ScreenDollarIndices
,
    D[f[t, r], r] + D[h[t, r], t] + PD[-Mx`a][v0[Mx`a]]
];

VerificationTest[
    0 == GCTensorCovDDiv[v1, 1, -decomp, chris][]
    - PDGChart[-decomp][-Mf`A][v1[Mf`A]]
    - chris[Mf`B, -Mf`B, -Mf`A] v1[Mf`A]
    // ContractGCTensors[#, holder1] & // NoScalar // ToCanonical[#, UseMetricOnVBundle -> None] &
];

VerificationTest[
    GCTensorCovDDiv[t1, 1, -decomp, chris] - ETensor[ContractGCTensors[
        PDGChart[-decomp][-Mf`B][t1[Mf`B, -Mf`A]]
        + chris[Mf`B, -Mf`B, -Mf`C] t1[Mf`C, -Mf`A]
        - chris[Mf`C, -Mf`B, -Mf`A] t1[Mf`B, -Mf`C]
    ,
        holder1
    ], {-Mf`A}] // NoScalar // ToCanonical[#, UseMetricOnVBundle -> None] & // ZeroGCTensorQ
];

VerificationTest[
    GCTensorCovDGrad[t1, -decomp, chris] - ETensor[ContractGCTensors[
        PDGChart[-decomp][-Mf`C][t1[Mf`A, -Mf`B]]
        + chris[Mf`A, -Mf`C, -Mf`D] t1[Mf`D, -Mf`B]
        - chris[Mf`D, -Mf`C, -Mf`B] t1[Mf`A, -Mf`D]
    ,
        holder1
    ], {Mf`A, -Mf`B, -Mf`C}] // NoScalar // ToCanonical[#, UseMetricOnVBundle -> None] & // ZeroGCTensorQ
];

VerificationTest[
    With[{
        term1 = ETensor[
            chris[Mf`D, -Mf`B, -Mf`E] p0[Mf`A, Mf`B, Mf`C, Mf`E]
            + chris[Mf`C, -Mf`B, -Mf`E]p0[Mf`A, Mf`B, Mf`E, Mf`D]
            + chris[Mf`B, -Mf`B, -Mf`E]p0[Mf`A, Mf`E, Mf`C, Mf`D]
            + chris[Mf`A, -Mf`B, -Mf`E]p0[Mf`E, Mf`B, Mf`C, Mf`D]
            + pd1[-Mf`B][p0[Mf`A, Mf`B, Mf`C, Mf`D]] // ContractGCTensors[holder1]
        ,
            {Mf`A, Mf`C, Mf`D}
        ],
        term2 = GCTensorCovDDiv[p0, 2, -decomp, chris]
    }, ToCanonical[term1 - term2, UseMetricOnVBundle -> None] // Simplify // ZeroGCTensorQ]
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["DeltaGCTensor"];

VerificationTest[
    With[{
        delta01 = DeltaGCTensor[decomp, -decomp],
        delta02 = GCTensorFromSparse[{{1, 1} -> 1, {2, 2} -> 1, {3, 3} -> ETensor[delta[Mx`a, -Mx`b]]}, {decomp, -decomp}]
    }, ZeroGCTensorQ[delta01 - delta02 // ToCanonical]]
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["Misc"];

VerificationTest[
    ContractGCTensors[n1[-Mf`A] n1[Mf`A], holder1, OtherReplaces -> { n1 -> GCTensorFromSparse[{{2} -> Sqrt[f[r]]}, {decomp}] }]
,
    1
];

MUnit`EndTestSection[];

MUnit`BeginTestSection["Change GCTensor indices"];

VerificationTest[
    With[{
        t1 = ContractGCTensors[CDMf[-Mf`A]@CDMf[Mf`A]@phi0[], holder1] /. PD[_]@phi0[] -> 0 // NoScalar // ToCanonical[#, UseMetricOnVBundle -> None] & // Simplify,
        t2 = 1 / Sqrt[detg] (-ParamD[t][Sqrt[detg] 1/h[r] ParamD[t]@phi0[]] + ParamD[r][Sqrt[detg] f[r] ParamD[r]@phi0[]])
    },
        Simplify[t1 - t2] == 0
    ]
];

GCTensor[{ETensor[n0[a], {a}]}, {NoDecomp, decomp}]

MUnit`EndTestSection[];

UndefTensor /@ {n0, n1, v0};
Undef /@ VisitorsOf@metricMf;
Undef /@ VisitorsOf@metricMx;
Undef /@ VisitorsOf@metricMr;
UndefMetric /@ {metricMx, metricMf, metricMr};
Undef /@ VisitorsOf@Mx;
Undef /@ VisitorsOf@Mf;
Undef /@ VisitorsOf@Mr;
UndefManifold /@ {Mx, Mf, Mr};
UndefConstantSymbol[dimx];
UndefConstantSymbol[L0];
UndefParameter /@ {r, t};
UndefScalarFunction /@ {h, f, Pn};
