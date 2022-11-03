BeginPackage["xTools`xDecomp`", {"xAct`xCore`", "xAct`xTensor`", "xTools`xTension`"}];

(Unprotect[#]; Remove[#];) & /@ Names@{"xTools`xDecomp`*", "xTools`xDecomp`Private`*"};

GChartQ::usage = "GChartQ[chart] gives True if chart is a generalized chart.";
CoordsOfGChart::usage = "CoordsOfGChart[chart] gives single-variable coordinate list of the generalized chart, each element is in the form {e, ed, param}."
SubManifoldsOfGChart::usage = "SubManifoldsOfGChart[chart] gives a list of basis tensors of the generalized chart.";
MetricOfGChart::usage = "MetricOfGChart[chart] gives the metric of generalized chart, if defined."
GChartMetricDecompose::usage = "GChartMetricDecompose[chart, metric[-a, -b]] decomposes the metric using the defined metric decomposition rule.";
VBundleOfGChart::usage = "VBundleOfGChart[chart] gives the base VBundle of the generalized chart.";
SubvbundlesOfGChart::usage = "SubvbundlesOfGChart[chart] gives the subvbundles of sub manifolds of the generalized chart.";
IndexDimensionOfGChart::usage = "IndexDimensionOfGChart[chart] gives the index dimension, i.e., the dimension where sub manifolds are treated as dim 1, of the chart.";

DefGBasis::usage = "DefGBasis[chart, {{e, ed, param}...}, {eT...}] defines a generalized chart with coodinate basis and sub manifolds."
SetGChartMetric::usage = "SetGChartMetric[chart, metric, metricValue, metricInvValue] defines metric decomposition rule for metric.";
DefGCovD::usage = "DefGCovD[covd, chart, metric, inverseMetric] defines a covd with the specified metric.";

GCovDQ::usage = "GCovDQ[covd] gives True for GCovD.";
GChartOfGCovD::usage = "GChartOfGCovD[covd] gives the GChart of the GCovD.";
MetricOfGCovD::usage = "MetricOfGCovD[covd] gives the metric defined for GCovD covd.";
InvMetricOfGCovD::usage = "InvMetricOfGCovD[covd] gives the inverse metric defined for GCovD covd.";

MetricGCovD::usage = "MetricGCovD[covd] represents the head of this metric defined to covd.";
PDGChart::usage = "PDGChart[chart, -a] represents the PD operator of the GChart chart.";
ExpandPDToBasis::usage = "ExpandPDToBasis[expr, chart] or ExpandPDToBasis[chart][expr] expands all partial derivatives in expr in terms of basis the coordinates and sub manifolds of chart.";
ExpandMetric::usage = "ExpandMetric[expr, chart] or ExpandMetric[chart][expr] replaces all metric in expr with the decomposed one defined in chart."
ExpandPDToGCTensor::usage = "ExpandPDToGCTensor[expr, chart] or ExpandPDToGCTensor[chart][expr] replaces all PDs to GCTensors in expr.";

GCTensor::usage = "GCTensor[values, charts] represents a generalized CTensor.";
GCTensorTranspose::usage = "GCTensorTranspose[GCTensor[...], perms] performs the tensor transpose of GCTensor.";
GCTensorContract::usage = "GCTensorContract[T, n1, n2] contracts the n1 axis with n2 axis of T.";
GCTensorContractTwo::usage = "GCTensorContractTwo[T1, T2, n1, n2] computes the tensor contraction between two tensors, in which n1-th axis of T1 and n2-th axis of T2.";
ContractGCTensors::usage = "ContractGCTensors[expr, covd] performs contractions of GCTensor's, use metric if needed.";
FixGCTensor::usage = "FixGCTensor[GCTensor[...]] fixes some ";
GCTensorPD::usage = "GCTensorPD[expr, chart] calculates the partial derivative of the tensor expression in chart.";

Begin["`Private`"];

GChartQ[c_] := CoordsOfGChart[c] =!= None;
CoordsOfGChart[_] = None;
CoordsOfGChart[-a_Symbol] := CoordsOfGChart[a];
SubManifoldsOfGChart[_] = None;
MetricOfGChart[_] = None;
GChartMetricDecompose[_, expr_] := expr;
VBundleOfGChart[chart_] := First@SlotsOfTensor@CoordsOfGChart[chart][[1, 1]];
SubvbundlesOfGChart[chart_] := SlotsOfTensor[#][[1]] & /@ SubManifoldsOfGChart[chart];
SyntaxInformation[GChartQ] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[CoordsOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SubManifoldsOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[MetricOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[GChartMetricDecompose] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[VBundleOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SubvbundlesOfGChart] = {"ArgumentsPattern" -> {_}};
Protect[GChartQ, CoordsOfGChart, CoordsOfGChart, SubManifoldsOfGChart, MetricOfGChart, GChartMetricDecompose, VBundleOfGChart, SubvbundlesOfGChart];

IndexDimensionOfGChart[chart_?GChartQ] := Length[CoordsOfGChart[chart]] + Length[SubManifoldsOfGChart[chart]];
SyntaxInformation[IndexDimensionOfGChart] = {"ArgumentsPattern" -> {_}};
Protect[IndexDimensionOfGChart];

DefGBasis[chart_Symbol, coordBasis_List, subManifolds_List] := Module[
    {params},
    GChartQ[chart] ^= True;
    CoordsOfGChart[chart] ^= {#[[1, 0]], #[[2, 0]], #[[3]]} & /@ coordBasis;
    SubManifoldsOfGChart[chart] ^= #[[0]] & /@ subManifolds;
    params = #[[3]] & /@ coordBasis;
    Replace[{e_[i1_Symbol], ed_[-i1_Symbol], param_} :> With[{
        M = BaseOfVBundle@VBundleOfIndex@i1,
        Q = Symbol[ToString[VBundleOfIndex@i1] <> "`Q"]
    },
        If[!xTensorQ[e],
            DefTensor[e[i1], M, PrintAs -> "(\!\(\*SubscriptBox[\(\[PartialD]\), \(" <> PrintAs[param] <> "\)]\))"]
        ];
        If[!xTensorQ[ed], DefTensor[ed[-i1], M, PrintAs -> "(d" <> PrintAs[param] <> ")"]];
        param /: PD[-a_Symbol?Q][param] := ed[-a];
        e /: e[a_Symbol] ed[-a_Symbol] = 1;
        e /: PD[_][e[a_Symbol]] = 0;
        e /: e[a_Symbol] PD[-a_Symbol][A_] := ParamD[param][A];
        ed /: PD[_][ed[-a_Symbol]] = 0;
    ]] /@ coordBasis;
    Replace[T_[i1_Symbol, -i2_Symbol] :> With[{
        M1 = BaseOfVBundle@VBundleOfIndex@i1,
        M2 = BaseOfVBundle@VBundleOfIndex@i2,
        Q = Symbol[ToString[VBundleOfIndex@i1] <> "`Q"]
    },
        If[!xTensorQ[T], DefTensor[T[i1, -i2], {M1, M2}]];
        T /: PD[_][T[-a_Symbol, b_Symbol]] = 0;
        T /: PD[_][T[a_Symbol, -b_Symbol]] = 0;
        T /: T[a_, b_Symbol] T[c_, -b_Symbol] := delta[a, c];
        T /: T[-a_Symbol, b_Symbol] PD[-b_Symbol][A_] := PD[-a][A];
        (# /: PD[-a_Symbol?Q]@# = 0) & /@ params;
    ]] /@ subManifolds;
    MapIndexed[With[{e1 = #1[[1, 0]], ed1 = #[[2, 0]], idx = #2[[1]]},
        With[{e2 = #[[1, 0]], ed2 = #[[2, 0]]},
            e1 /: e1[a_Symbol] ed2[-a_Symbol] = 0;
            e2 /: e2[a_Symbol] ed1[-a_Symbol] = 0;
        ] & /@ coordBasis[[idx + 1 ;;]];
    ] &, coordBasis];
    Outer[With[{e = #1[[1, 0]], ed = #1[[2, 0]], t = #2[[0]]},
        e /: e[a_Symbol] t[_, -a_Symbol] = 0;
        ed /: ed[-a_Symbol] t[_, a_Symbol] = 0;
    ] &, coordBasis, subManifolds, 1];
    MapIndexed[With[{t1 = #1[[0]], idx = #2[[1]]},
    With[{t2 = #[[0]]},
        t1 /: t1[_, a_Symbol] t2[_, -a_Symbol] = 0;
        t1 /: t1[_, -a_Symbol] t2[_, a_Symbol] = 0;
        ] & /@ subManifolds[[idx + 1 ;;]];
    ] &, subManifolds];
];
SyntaxInformation[DefGBasis] = {"ArgumentsPattern" -> {_, _, _}};
Protect[DefGBasis];

Quiet[
    IndexToPattern[idx_Symbol] := Pattern[idx, Blank[]];
    IndexToPattern[-idx_Symbol] := -Pattern[idx, Blank[]],
{RuleDelayed::rhs}];
SetGChartMetric[chart_, metric_, mval_, inv_] := With[{
    metricPat = metric @@ IndexToPattern /@ FindFreeIndices[mval],
    metricInvPat = metric @@ IndexToPattern /@ FindFreeIndices[inv]
},
    MetricOfGChart[chart] ^= metric;
    With[{inds = List @@ DeleteDuplicates[UpIndex /@ FindDummyIndices@mval]},
        chart /: GChartMetricDecompose[chart, metricPat] := Module[inds, mval]
    ];
    With[{inds = List @@ DeleteDuplicates[UpIndex /@ FindDummyIndices@inv]},
        chart /: GChartMetricDecompose[chart, metricInvPat] := Module[inds, inv]
    ];
];
SyntaxInformation[SetGChartMetric] = {"ArgumentsPattern" -> {_, _, _, _}};
Protect[SetGChartMetric];

GCovDQ[e_] := GChartOfGCovD[e] =!= None;
GChartOfGCovD[_] = None;
MetricOfGCovD[_] = None;
InvMetricOfGCovD[_] = None;
SyntaxInformation[GCovDQ] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[GChartOfGCovD] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[MetricOfGCovD] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[InvMetricOfGCovD] = {"ArgumentsPattern" -> {_}};
Protect[GCovDQ, GChartOfGCovD, MetricOfGCovD, InvMetricOfGCovD];

DefGCovD[covd_Symbol, chart_?GChartQ, metric_, metricInv_] := Module[
    {},
    GChartOfGCovD[covd] ^= chart;
    MetricOfGCovD[covd] ^= metric;
    InvMetricOfGCovD[covd] ^= metricInv;
];
SyntaxInformation[DefGCovD] = {"ArgumentsPattern" -> {_, _, _, _}};
Protect[DefGCovD];

MetricGCovD[covd_?GCovDQ][-a_Symbol, -b_Symbol] := MetricOfGCovD[covd][-a, -b];
MetricGCovD[covd_?GCovDQ][a_Symbol, b_Symbol] := InvMetricOfGCovD[covd][a, b];
SyntaxInformation[MetricGCovD] = {"ArgumentsPattern" -> {_}};
Protect[MetricGCovD];

PDGChart[chart_?GChartQ, -a_Symbol][t_GCTensor[inds__]] := GCTensorPD[t, chart][-a, inds];
SyntaxInformation[PDGChart] = {"ArgumentsPattern" -> {_, _}};
Protect[PDGChart];

ExpandPDToBasis[chart_?GChartQ][expr_] := ExpandPDToBasis[expr, chart];
ExpandPDToBasis[expr_, chart_?GChartQ] := Module[
    {vb, a, A, replacement1, replacement2, rule},
    vb = VBundleOfGChart@chart;
    a = GetIndicesOfVBundle[vb, 1][[1]];
    replacement1[a_, A_] = Plus @@ (ParamD[#[[3]]][A] #[[2]][-a] & /@ CoordsOfGChart[chart]);
    replacement2[a_, A_] = With[{
        i1 = GetIndicesOfVBundle[First@SlotsOfTensor@#, 1][[1]],
        e1 = #
    }, {{i1}, e1[i1, -a] PD[-i1][A]}] & /@ SubManifoldsOfGChart[chart];
    rule = With[{
        Q = Symbol[ToString[vb] <> "`Q"]
    }, PD[-a_Symbol?Q][A_] :> replacement1[a, A] + Plus @@ Module @@@ replacement2[a, A]];
    expr /. rule
];
SyntaxInformation[ExpandPDToBasis] = {"ArgumentsPattern" -> {_, _.}};
Protect[ExpandPDToBasis];

ExpandPDToGCTensor[chart_][expr_] := ExpandPDToGCTensor[expr, chart];
ExpandPDToGCTensor[expr_, chart_?GChartQ] := With[{
    Q = Symbol[ToString[VBundleOfGChart@chart] <> "`Q"]
}, expr /. PD[-a_Symbol?Q] -> PDGChart[chart, -a]];
SyntaxInformation[ExpandPDToGCTensor] = {"ArgumentsPattern" -> {_, _.}};
Protect[ExpandPDToGCTensor];

ExpandMetric[chart_?GChartQ][expr_] := ExpandMetric[expr, chart];
ExpandMetric[expr_, chart_?GChartQ] := With[{metric = MetricOfGChart@chart},
   expr /. metric[inds__] :> GChartMetricDecompose[chart, metric[inds]]
];
SyntaxInformation[ExpandMetric] = {"ArgumentsPattern" -> {_, _.}};
Protect[ExpandMetric];

(* GCTensor *)
MoveTo[n_, n_] := {1};
MoveTo[from_Integer, to_Integer] := If[
    from < to,
    Join[Range[from - 1], {to}, Range[from, to - 1]],
    Join[Range[to - 1], Range[to + 1, from], {to}]
];
UpBasis[a_Symbol] := a;
UpBasis[-a_Symbol] := a;
CompletePerm[perm_, len_] := With[{
    l = Length@perm
}, If[len > l, Join[perm, Range[l + 1, len]], perm]];

GCTensorDataTranspose[array_, basis_, perms_] := With[{
    subMStartI = Length[CoordsOfGChart@UpBasis@#] & /@ basis,
    cperms = CompletePerm[perms, Length@basis]
},
    MapIndexed[Function[{elem, inds},
        ETensorTranspose[elem, InversePermutation@Ordering[
            cperms[[
                #[[1]] & /@ Position[inds - subMStartI, _?(# > 0 &), {1}]
            ]]
        ]]
    ], array, {Length@basis}] // Transpose[#, perms] &
];

GCTensorTranspose[GCTensor[array_, basis_], perms_] := (
    GCTensorDataTranspose[array, basis, perms] // GCTensor[#, Permute[basis, perms]] &
);
SyntaxInformation[GCTensorTranspose] = {"ArgumentsPattern" -> {_, _}};
Protect[GCTensorTranspose];

DualIndexQ[-a_Symbol, a_Symbol] = True;
DualIndexQ[a_Symbol, -a_Symbol] = True;
DualIndexQ[_, _] = False;
ContractableChartsQ[c_Symobl, -c_Symbol] = True;
ContractableChartsQ[-c_Symobl, c_Symbol] = True;
ContractableChartsQ[-c1_Symbol, c2_Symbol] := ContractableChartsQ[c1, -c2];
ContractableChartsQ[c1_Symbol, -c2_Symbol] := Length[CoordsOfGChart[c1]] === Length[CoordsOfGChart[c2]] && Length[SubManifoldsOfGChart[c1]] === Length[SubManifoldsOfGChart[c2]];
ContractableChartsQ[_, _] = False;

Null@ContractPair;

ContractKernels[basis1_, basis2_] := Module[
    {ai, leftI, rightI, leftAi, rightAi, cbasisInners, subMInners},
    ai = UniqueIndex@First@GetIndicesOfVBundle[VBundleOfGChart@UpBasis@basis1, 1];
    {leftI, rightI} = If[MatchQ[basis2, _Times], {1, 2}, {2, 1}];
    leftAi = (3 - 2*leftI) * ai;
    rightAi = (3 - 2*rightI) * ai;
    cbasisInners = #1[[leftI]][leftAi] #2[[rightI]][rightAi] & @@@ Thread@{CoordsOfGChart@UpBasis@basis1, CoordsOfGChart@UpBasis@basis2};
    subMInners = Module[
        {leftA, rightA},
        leftA = (2*leftI - 3) * UpIndex[First@GetIndicesOfVBundle[SlotsOfTensor[#1][[1]], 1]];
        rightA = (2*rightI - 3) * UpIndex[First@GetIndicesOfVBundle[SlotsOfTensor[#2][[1]], 1, {UpIndex@leftA}]];
        ETensor[#1[leftA, leftAi] #2[rightA, rightAi], {leftA, rightA}]
    ] & @@@ Thread@{SubManifoldsOfGChart@UpBasis@basis1, SubManifoldsOfGChart@UpBasis@basis2};
    {cbasisInners, subMInners}
]

ApplyContractKernel[kernels_, expr_] := With[{
    tensors = Cases[kernels, _ETensor],
    others = Cases[kernels, a_ /; !MatchQ[a, _ETensor]]
},
    Fold[ETensorContractTwo[#2, #1, {2}, {Length@tensors}] &, expr, Reverse@tensors] * (Times @@ others)
];
GCTensorDotInner[lhs_, ContractPair[rhs_, len_]] := With[{r1 = Range[-len, -1], r2 = Range[len]},
    ETensorContractTwo[lhs, rhs, r1, r2]
];
GCTensorDot[GCTensor[arr1_, basis1_], GCTensor[arr2_, basis2_], len_] := Module[
    {ckernels, coordLengths, a1, a2},
    ckernels = (Join @@ ContractKernels[#1, #2]) & @@@ Thread@{basis1[[Length[basis1] - len + 1 ;;]], Reverse@basis2[[;; len]]};
    coordLengths = Length@CoordsOfGChart@UpBasis@# & /@ basis1[[Length[basis1] - len + 1 ;;]];
    a1 = If[len > 1, Map[Flatten, arr1, {Length[basis1] - len}], arr1];
    a2 = MapIndexed[Function[{elem, indices}, With[{
            selectedKernels = #1[[#2]] & @@@ Thread@{ckernels, indices},
            contractLen = Length@Cases[Thread[indices > coordLengths], True]
        }, Map[ContractPair[
            ApplyContractKernel[selectedKernels, #],
            contractLen
        ] &, elem, {Length@basis2 - len}]]
    ], arr2, {len}] // If[len > 1, Flatten[#, len - 1], #] &;
    GCTensor[
        Inner[GCTensorDotInner, a1, a2, Plus],
        Join[
            basis1[[;; Length[basis1] - len]],
            basis2[[len + 1 ;;]]
        ]
    ]
];

GCTensorContractTwo::icpbs = "cannot contract incompatible basis `1` and `2`.";
GCTensorContractTwo[GCTensor[arr1_, basis1_], GCTensor[arr2_, basis2_], n1_List, n2_List] := Module[
    {d1, d2, perm1, perm2, trans1, trans2},
    d1 = Length@basis1;
    d2 = Length@basis2;
    If[!ContractableChartsQ[#1, #2], Throw@Message[GCTensorContractTwo::icpbs, #1, #2]] & @@@ Thread@{basis1[[#]] & /@ n1, basis2[[#]] & /@ n2};
    perm1 = PermutationProduct @@ MapIndexed[CompletePerm[MoveTo[#1, d1 - #2[[1]] + 1], d1] &, Reverse@n1];
    perm2 = PermutationProduct @@ MapIndexed[CompletePerm[MoveTo[#1, #2[[1]]], d2] &, n2];
    trans1 = GCTensorTranspose[GCTensor[arr1, basis1], perm1];
    trans2 = GCTensorTranspose[GCTensor[arr2, basis2], perm2];
    GCTensorDot[trans1, trans2, Length@n1]
];
SyntaxInformation[GCTensorContractTwo] = {"ArgumentsPattern" -> {_, _, _, _}};
Protect[GCTensorContractTwo];

GCTensorContractLast2[GCTensor[arr_, basis_]] := Module[
    {cbasis, csubMs, ccount},
    ccount = Length@CoordsOfGChart@UpBasis@basis[[-1]];
    {cbasis, csubMs} = ContractKernels[basis[[-1]], basis[[-2]]];
    Map[Function[elem, With[{diag = Diagonal@elem},
        Total[diag[[;; ccount]] * cbasis] + Total[ETensorContractTwo[#1, #2, {-2, -1}, {1, 2}] & @@@ Thread@{diag[[ccount + 1 ;;]], csubMs}]
    ]], arr, {Length@basis - 2}] // GCTensor[#, basis[[;; Length@basis - 2]]] &
];

GCTensorContract::icpbs = "cannot contract incompatible basis `1` and `2`.";
GCTensorContract[t_, n1_List, n2_List] := Fold[GCTensorContract[#1, #2[[1]], #2[[2]]], t, Thread@{n1, n2}];
GCTensorContract[GCTensor[arr_, basis_], n1_Integer, n2_Integer] := Module[
    {d, perm},
    d = Length@basis;
    If[!ContractableChartsQ[#1, #2], Throw@Message[GCTensorContract::icpbs, #1, #2]] &[basis[[n1]], basis[[n2]]];
    perm = PermutationProduct[CompletePerm[MoveTo[n1, d], d], CompletePerm[MoveTo[n2, d - 1], d]];
    GCTensorContractLast2[GCTensorTranspose[GCTensor[arr, basis], perm]]
];
SyntaxInformation[GCTensorContract] = {"ArgumentsPattern" -> {_, _, _}};
Protect[GCTensorContract];

DefGCTensorMapFunc[funcs__] := Function[func,
    GCTensor /: func[GCTensor[arr_, basis_], args___] := GCTensor[Map[func[#, args] &, arr, {Length@basis}], basis];
    GCTensor /: func[GCTensor[arr_, basis_][inds__], args___] := GCTensor[Map[func[#, args] &, arr, {Length@basis}], basis][inds];
] /@ {funcs};

DefGCTensorMapFunc[
    ScreenDollarIndices,
    ToCanonical,
    ChangeCovD,
    ChangeCurvature,
    ChristoffelToRiemann,
    ChristoffelToGradMetric
];
GCTensor[e_, {}][] := Scalar[e];
GCTensor /: ParamD[args__]@GCTensor[arr_, basis] := GCTensor[Map[ParamD[args], arr, {Length@basis}], basis];
GCTensor /: ETensor[t_GCTensor[inds__], {inds2__}] := With[{
    perm = Ordering[{inds2}][[Ordering@Ordering@{inds}]]
}, GCTensorTranspose[t, perm]] /; Sort@{inds} === Sort@{inds2};

GCTensor /: GCTensor[arr1_, basis_] + GCTensor[arr2_, basis_] := GCTensor[arr1 + arr2, basis];
GCTensor /: x_?xTools`xTension`Private`IndexedScalarQ * GCTensor[arr_, basis_] := GCTensor[arr * x, basis];
GCTensor /: GCTensor[arr1_, basis1_][inds1__] + GCTensor[arr2_, basis2_][inds2__] := With[{
    perm = Ordering[{inds1}][[Ordering@Ordering@{inds2}]]
},
    GCTensor[arr1 + GCTensorDataTranspose[arr2, basis2, perm], basis1][inds1]
] /; Sort@{inds1} === Sort@{inds2} && Length@basis1 === Length@basis2;

GCTensor /: x_?xAct`xTensor`Private`NonIndexedScalarQ * GCTensor[arr_, basis_][inds__] := GCTensor[arr * x, basis][inds];

SyntaxInformation[GCTensor] = {"ArgumentsPattern" -> {_, _}};
Protect[GCTensor];

GCTensorPD[chart_][expr_] := GCTensorPD[expr, chart];
GCTensorPD[GCTensor[arr_, basis_], chart_] := With[{
    paramds = ParamD[#[[3]]] & /@ CoordsOfGChart@chart,
    subvbs = ETensorPD@First@SlotsOfTensor@# & /@ SubManifoldsOfGChart@chart,
    len = Length@basis
},
    GCTensor[Join[
        Map[#, arr, {len}] & /@ paramds,
        Map[#, arr, {len}] & /@ subvbs
    ], Prepend[basis, -chart]]
];
SyntaxInformation[GCTensorPD] = {"ArgumentsPattern" -> {_, _.}};
Protect[GCTensorPD];

ContractTwoIndexedGCTensors[t1_GCTensor[inds1__], t2_GCTensor[inds2__]] := Module[
    {pairs, pos1, pos2, res},
    pairs = xAct`xTensor`Private`TakePairs[{inds1}, {inds2}];
    pos1 = FirstPosition[{inds1}, #][[1]] & /@ pairs;
    pos2 = FirstPosition[{inds2}, ChangeIndex@#][[1]] & /@ pairs;
    res = GCTensorContractTwo[t1, t2, pos1, pos2];
    res @@ Join[Delete[{inds1}, Transpose@{pos1}], Delete[{inds2}, Transpose@{pos2}]]
];

ContractOneIndexedGCTensors[t_GCTensor[inds__]] := Module[
    {pairs, pos1, pos2, res},
    pairs = xAct`xTensor`Private`TakePairs@{inds};
    pos1 = FirstPosition[{inds}, #][[1]] & /@ pairs;
    pos2 = FirstPosition[{inds}, ChangeIndex@#][[1]] & /@ pairs;
    res = GCTensorContract[t, pos1, pos2];
    res @@ Delete[{inds}, Join[Transpose@{pos1}, Transpose@{pos2}]]
];

ExecuteContraction[{_, idx_, idx_}, l_List] := {
    ContractOneIndexedGCTensors[l[[idx]]],
    Delete[l, idx]
};
ExecuteContraction[{_, idx1_, idx2_}, l_List] := {
    ContractTwoIndexedGCTensors[l[[idx1]], l[[idx2]]],
    Delete[l, {{idx1}, {idx2}}]
} /; idx1 != idx2;
OptimizedGCTensorContraction[{}] = 1;
OptimizedGCTensorContraction[l_List] := Module[
    {tensors, contractions, contractionCounts, res, l2},
    tensors = MapIndexed[{List @@ #1, #2[[1]]} &, l];
    contractions = Flatten[
        Outer[{Length@xAct`xTensor`Private`DropPairs[#1[[1]], #2[[1]]], #1[[2]], #2[[2]]} &, tensors, tensors, 1],
        1
    ];
    contractionCounts = #[[1]] & /@ contractions;
    {res, l2} = ExecuteContraction[contractions[[FirstPosition[contractionCounts, Min@contractionCounts][[1]]]], l];
    res * OptimizedGCTensorContraction[l2]
];

ContractGCTensors[covd_][expr_] := ContractGCTensors[expr, covd];
ContractGCTensors[expr_Plus, covd_] := ContractGCTensors[#, covd] & /@ expr;
ContractGCTensors[expr_Times, covd_] := With[{
    tensors = Cases[expr, GCTensor[__][__]],
    others = ContractGCTensors[#, covd] & /@ DeleteCases[List @@ expr, GCTensor[__][__]]
}, (Times @@ others) * OptimizedGCTensorContraction[tensors]];
ContractGCTensors[expr_, _] := expr;
SyntaxInformation[ContractGCTensors] = {"ArgumentsPattern" -> {_, _.}};
Protect[ContractGCTensors];

End[];

EndPackage[];