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

ExpandPD::usage = "ExpandPD[expr, chart] or ExpandPD[chart][expr] expands all partial derivatives in expr into the coordinates and sub manifolds of chart.";
ExpandMetric::usage = "ExpandMetric[expr, chart] or ExpandMetric[chart][expr] replaces all metric in expr with the decomposed one defined in chart."

GCTensor::usage = "GCTensor[values, charts] represents a generalized CTensor.";
GCTensorTranspose::usage = "GCTensorTranspose[GCTensor[...], perms] performs the tensor transpose of GCTensor.";
GCTensorContractTwo::usage = "GCTensorContractTwo[T1, T2, n1, n2] computes the tensor product and contracts n1-th axis of T1 and n2-th axis of T2.";
FixGCTensor::usage = "FixGCTensor[GCTensor[...]] fixes some ";

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
    {},
    GChartQ[chart] ^= True;
    CoordsOfGChart[chart] ^= {#[[1, 0]], #[[2, 0]], #[[3]]} & /@ coordBasis;
    SubManifoldsOfGChart[chart] ^= #[[0]] & /@ subManifolds;
    Replace[{e_[i1_Symbol], ed_[-i1_Symbol], param_} :> With[{
        M = BaseOfVBundle@VBundleOfIndex@i1
    },
        If[!xTensorQ[e],
            DefTensor[e[i1], M, PrintAs -> "(\!\(\*SubscriptBox[\(\[PartialD]\), \(" <> PrintAs[param] <> "\)]\))"]
        ];
        If[!xTensorQ[ed], DefTensor[ed[-i1], M, PrintAs -> "(d" <> PrintAs[param] <> ")"]];
        param /: PD[-a_Symbol][param] := ed[-a];
        e /: e[a_Symbol] ed[-a_Symbol] = 1;
        e /: PD[_][e[a_Symbol]] = 0;
        e /: e[a_Symbol] PD[-a_Symbol][A_] := ParamD[param][A];
        ed /: PD[_][ed[-a_Symbol]] = 0;
    ]] /@ coordBasis;
    Replace[T_[i1_Symbol, -i2_Symbol] :> With[{
        M1 = BaseOfVBundle@VBundleOfIndex@i1,
        M2 = BaseOfVBundle@VBundleOfIndex@i2
    },
        If[!xTensorQ[T], DefTensor[T[i1, -i2], {M1, M2}]];
        T /: PD[_][T[-a_Symbol, b_Symbol]] = 0;
        T /: PD[_][T[a_Symbol, -b_Symbol]] = 0;
        T /: T[a_, b_Symbol] T[c_, -b_Symbol] := delta[a, c];
        T /: T[-a_Symbol, b_Symbol] PD[-b_Symbol][A_] := PD[-a][A];
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

ExpandPD[chart_?GChartQ][expr_] := ExpandPD[expr, chart];
ExpandPD[expr_, chart_?GChartQ] := Module[
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
SyntaxInformation[ExpandPD] = {"ArgumentsPattern" -> {_, _.}};
Protect[ExpandPD];

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

GCTensorTranspose[GCTensor[array_, basis_], perms_] := With[{
    subMStartI = Length[CoordsOfGChart@UpBasis@#] & /@ basis,
    cperms = CompletePerm[perms, Length@basis]
},
    MapIndexed[Function[{elem, inds},
        ETensorTranspose[elem, Ordering[
            cperms[[
                #[[1]] & /@ Position[inds - subMStartI, _?(# > 0 &), {1}]
            ]]
        ]]
    ], array, {Length@basis}] // GCTensor[Transpose[#, perms], Permute[basis, perms]] &
];
SyntaxInformation[GCTensorTranspose] = {"ArgumentsPattern" -> {_, _}};
Protect[GCTensorTranspose];

DualIndexQ[-a_Symbol, a_Symbol] = True;
DualIndexQ[a_Symbol, -a_Symbol] = True;
DualIndexQ[_, _] = False;
ContractableChartsQ[c_Symobl, -c_Symbol] = True;
ContractableChartsQ[-c_Symobl, c_Symbol] = True;
ContractableChartsQ[-c1_Symbol, c2_Symbol] := ContractableChartsQ[c1, -c2];
ContractableChartsQ[c1_Symbol, -c2_Symbol] := Length[CoordsOfGChart[c1]] === Length[CoordsOfGChart[c2]] && Length[SubManifoldsOfGChart[c1]] === Length[SubManifoldsOfGChart[c2]];

Null@ContractPair;

GCTensorDot[GCTensor[arr1_, basis1_], GCTensor[arr2_, basis2_]] := Module[
    {b1, b2, leftI, rightI, ai, leftAi, rightAi, cbasisInners, subMInners},
    b1 = basis1[[-1]];
    b2 = basis2[[1]];
    ai = UniqueIndex@First@GetIndicesOfVBundle[VBundleOfGChart@UpBasis@b1, 1];
    {leftI, rightI} = If[MatchQ[b2, _Times], {1, 2}, {2, 1}];
    leftAi = (3 - 2*leftI) * ai;
    rightAi = (3 - 2*rightI) * ai;
    cbasisInners = #1[[leftI]][leftAi] #2[[rightI]][rightAi] & @@@ Thread@{CoordsOfGChart@UpBasis@b1, CoordsOfGChart@UpBasis@b2};
    subMInners = Module[
        {leftA, rightA},
        leftA = (2*leftI - 3) * UpIndex[First@GetIndicesOfVBundle[SlotsOfTensor[#1][[1]], 1]];
        rightA = (2*rightI - 3) * UpIndex[First@GetIndicesOfVBundle[SlotsOfTensor[#2][[1]], 1, {leftA}]];
        ETensor[#1[leftA, leftAi] #2[rightA, rightAi], {leftA, rightA}]
    ] & @@@ Thread@{SubManifoldsOfGChart@UpBasis@b1, SubManifoldsOfGChart@UpBasis@b2};

    Inner[If[#2[[2]], ETensorDot[#1, #2[[1]]], ETensorProduct[#1, #2[[1]]]] &, arr1, Thread@ContractPair[
        Join[
            cbasisInners * arr2[[1 ;; Length@cbasisInners]],
            ETensorDot @@@ Thread@{subMInners, arr2[[Length[cbasisInners] + 1 ;;]]}
        ],
        Array[# > Length@cbasisInners &, Length[cbasisInners] + Length[subMInners]]
    ], Plus]
];

GCTensorContractTwo::icpbs = "cannot contract incompatible basis `1` and `2`.";
GCTensorContractTwo[GCTensor[arr1_, basis1_], GCTensor[arr2_, basis2_], n1_Integer, n2_Integer] := With[
    {b1 = basis1[[n1]], b2 = basis2[[n2]]},
    If[!ContractableChartsQ[b1, b2], Throw@Message[GCTensorContractTwo::icpbs, b1, b2]];
    GCTensorDot[
        GCTensorTranspose[GCTensor[arr1, basis1], MoveTo[n1, Length@basis1]],
        GCTensorTranspose[GCTensor[arr2, basis2], MoveTo[n2, 1]]
    ]
];
SyntaxInformation[GCTensorContractTwo] = {"ArgumentsPattern" -> {_, _, _, _}};
Protect[GCTensorContractTwo];

End[];

EndPackage[];