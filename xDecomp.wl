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
GetAllHeldTensors::usage = "GetAllCachedTensors[holder] gives a list of tensors that has a definition with the holder.";
GetAllHeldTensorRules::usage = "GetAllHeldTensorRules[holder] gives a rule list to replace all tensors defined in holder as the values.";
ClearGCTensorHolderCache::usage = "ClearGCTensorHolderCache[holder] clears the tensor cache.";

GCTensorHolderQ::usage = "GCTensorHolderQ[covd] gives True for GCovD.";

SetHeldMetric::usage = "SetHeldMetric[holder, metric, invMetric] sets the metric of the holder.";
AddCurvatureTensorsToHolder::usage = "AddCurvatureTensorsToHolder[holder, chart, christoffelTensor] defines curvature tensors to the holder, using the given Christoffel tensor.";

HeldGCTensor::usage = "HeldGCTensor[holder, tensor] gives the tensor provider of the tensor.";
GCTensorHolderAction::usage = "GCTensorHolderAction[holder, tag, expr] defines some extra action applies in various computation stages. Typically used to simplify expressions.";
GCTensorHolderDefaultAction::usage = "GCTensorHolderDefaultAction[holder, tag, expr] is the default value for GCTensorHolderAction[holder, tag, expr].";
CachedGCTensor::usage = "CachedGCTensor[holder, tensor, inds] gives the cached tensor with specific risen/lowered indices.";

PDGChart::usage = "PDGChart[chart, -a] represents the PD operator of the GChart chart.";
ExpandPDToBasis::usage = "ExpandPDToBasis[expr, chart] or ExpandPDToBasis[chart][expr] expands all partial derivatives in expr in terms of basis the coordinates and sub manifolds of chart.";
ExpandMetric::usage = "ExpandMetric[expr, chart] or ExpandMetric[chart][expr] replaces all metric in expr with the decomposed one defined in chart."
ExpandPDToGCTensor::usage = "ExpandPDToGCTensor[expr, chart] or ExpandPDToGCTensor[chart][expr] replaces all PDs to GCTensors in expr.";

GCTensor::usage = "GCTensor[values, charts] represents a generalized CTensor.";
GCTensorTranspose::usage = "GCTensorTranspose[GCTensor[...], perms] performs the tensor transpose of GCTensor.";
GCTensorContract::usage = "GCTensorContract[T, n1, n2] contracts the n1 axis with n2 axis of T.";
GCTensorContractTwo::usage = "GCTensorContractTwo[T1, T2, n1, n2] computes the tensor contraction between two tensors, in which n1-th axis of T1 and n2-th axis of T2.";
GCTensorProduct::usage = "GCTensorProduct[T1, T2] computes the tensor product of T1 and T2.";
GCTensorFixedContract::usage = "GCTensorFixedContract[T, T2, n] contracts n-th axis of T with the symmetric tensor T2, with the resulting axis stay at the original position. Typically used to change index.";
GCTensorToBasis::usage = "GCTensorToBasis[expr] converts GCTensor expressions to basis.";

GCTensorToComponentRules::usage = "GCTensorToComponentRules[T, head] converts GCTensor to a list of component rules.";
HeadOfCoordinateIndex::usage = "HeadOfCoordinateIndex is an option of GCTensorToComponentRules that defines the head used for coordinate indices. Default is LI.";
ExcludeZeros::usage = "ExcludeZeros is a boolean option of GCTensorToComponentRules. When True, only non-zero components are listed.";

ContractGCTensors::usage = "ContractGCTensors[expr, covd] performs contractions of GCTensor's, use metric if needed.";
ValidateGCTensor::usage = "ValidateGCTensor[GCTensor[...]] fixes some ";
ZeroGCTensor::usage = "ZeroGCTensor[charts] creats a GCTensor with zero components.";
CreateGCTensor::usage = "CreateGCTensor[{{r1, r2, ...} -> expr}, charts] is a convenient function for creating GCTensors.";
GCTensorPD::usage = "GCTensorPD[expr, chart] calculates the partial derivative of the tensor expression in chart.";

RecoverSubRiemannTensors::usage = "RecoverSubRiemannTensors[GCTensor[..., {-c1, -c2, -c3, -c4}]] replaces Christoffel tensors of sub manifolds with Riemann tensors.";

$xDecompVerbose = False;
$xDecompVerbose::usage = "$xDecompVerbose is a global boolean variable, when set to True, messages wil be printed.";

Begin["`Private`"];

Off[RuleDelayed::rhs];

GChartQ[c_] := CoordsOfGChart[c] =!= None;
CoordsOfGChart[_] = None;
CoordsOfGChart[-a_Symbol] := CoordsOfGChart[a];
SubManifoldsOfGChart[_] = None;
MetricOfGChart[_] = None;
GChartMetricDecompose[_, expr_] := expr;
VBundleOfGChart[chart_] := First@SlotsOfTensor@CoordsOfGChart[chart][[1, 1]];
SubvbundlesOfGChart[chart_Symbol] := SlotsOfTensor[#][[1]] & /@ SubManifoldsOfGChart[chart];
SubvbundlesOfGChart[-chart_Symbol] := -SubvbundlesOfGChart[chart];
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

GCTensorHolderQ[_] = None;
SyntaxInformation[GCTensorHolderQ] = {"ArgumentsPattern" -> {_}};
Protect[GCTensorHolderQ];

UpDownIndexToNumber[_Symbol] = 1;
UpDownIndexToNumber[-_Symbol] = -1;
UpDownIndexToNumber[_] = 0;

HeldGCTensor::undef = "tensor provider for `1` of holder `2` is not defined.";
HeldGCTensor[holder_, tensor_] := Throw@Message[HeldGCTensor::undef, tensor, holder];
SyntaxInformation[HeldGCTensor] = {"ArgumentsPattern" -> {_, _}};
Protect[HeldGCTensor];

GCTensorHolderAction[holder_, tag_, expr_] := GCTensorHolderDefaultAction[holder, tag, expr];
GCTensorHolderDefaultAction[holder, _, e_] := e;
GCTensorHolderDefaultAction[_, "PostChangeIndex"[tensor_, inds_, i_], expr_] := Simplification@ContractMetric[expr];
GCTensorHolderDefaultAction[_, "PostChangeIndex"[tensor_, inds_, i_], expr: GCTensor[arr_, basis_]] := Simplification@ContractMetric[expr, First@MetricsOfVBundle@# & /@ SubvbundlesOfGChart@UpBasis@basis[[i]]];
GCTensorHolderDefaultAction[_, "PostCurvatureTensorCalculation"[tensor_], expr_] := Simplify@ToCanonical[expr, UseMetricOnVBundle -> None];
SyntaxInformation[GCTensorHolderAction] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[GCTensorHolderDefaultAction] = {"ArgumentsPattern" -> {_, _, _}};
Protect[GCTensorHolderAction, GCTensorHolderDefaultAction];

CachedGCTensor[holder_Symbol, tensor_][inds___] := CachedGCTensor[holder, tensor, UpDownIndexToNumber /@ {inds}][inds];
CachedGCTensor[holder_Symbol, tensor_, inds_] := (
    If[$xDecompVerbose, Print["[CachedGCTensor] calculating tensor ", tensor]];
    holder /: CachedGCTensor[holder, tensor, inds] = HeldGCTensor[holder, tensor]
) /; inds === UpDownIndexToNumber /@ SlotsOfTensor@tensor;
CachedGCTensor[holder_Symbol, tensor_, inds_] := With[{
    firstUnmatched = First@FirstPosition[MapThread[#1 =!= UpDownIndexToNumber@#2 &, {inds, SlotsOfTensor@tensor}], True],
    slots = SlotsOfTensor@tensor
},
    holder /: CachedGCTensor[holder, tensor, inds] = With[{
        tensor2 = CachedGCTensor[holder, tensor, MapAt[-# &, inds, firstUnmatched]]
    },
        If[$xDecompVerbose, Print["[CachedGCTensor] ", "calculating index change ", tensor, inds, " from ", MapAt[-# &, inds, firstUnmatched]]];
        GCTensorFixedContract[
            tensor2,
            CachedGCTensor[holder, First@MetricsOfVBundle@slots[[firstUnmatched]], {#, #} &[inds[[firstUnmatched]]]],
            firstUnmatched
        ] // GCTensorHolderAction[holder, "PostChangeIndex"[tensor, inds, firstUnmatched], #] &
    ]
] /; inds =!= UpDownIndexToNumber /@ SlotsOfTensor@tensor;
SyntaxInformation[CachedGCTensor] = {"ArgumentsPattern" -> {_, _, _.}};
Protect[CachedGCTensor];

PersistentTensorCacheQ[_, _] = False;
SetHeldMetric[holder_Symbol, metric_, metrici_GCTensor, invMetric_GCTensor] := (
    holder /: CachedGCTensor[holder, metric, {-1, -1}] = metrici;
    holder /: CachedGCTensor[holder, metric, {1, 1}] = invMetric;
    holder /: PersistentTensorCacheQ[holder, metric] = True;
);
SyntaxInformation[SetHeldMetric] = {"ArgumentsPattern" -> {_, _, _, _}};
Protect[SetHeldMetric];

AddCurvatureTensorsToHolder[holder_Symbol, chart_?GChartQ, chris_] := Module[
    {vb, metric, cd, temp, a0, b0, c0, d0},
    vb = VBundleOfGChart@chart;
    metric = First@MetricsOfVBundle@vb;
    cd = MasterOf@chris;
    {a0, b0, c0, d0} = GetIndicesOfVBundle[vb, 4];

    With[{
        a = a0,
        b = b0,
        c = c0,
        d = d0,
        subCDs = CovDOfMetric@First@MetricsOfVBundle@# & /@ SubvbundlesOfGChart@chart,
        subMetrics = First@MetricsOfVBundle@# & /@ SubvbundlesOfGChart@chart,
        metric = First@MetricsOfVBundle@vb,
        cd2 = cd,
        riem = Riemann@cd,
        ricci = Ricci@cd,
        ricciScalar = RicciScalar@cd
    },
        (* Christoffel tensor *)
        holder /: HeldGCTensor[holder, chris] := ETensor[
            ChristoffelToGradMetric[chris[a, -b, -c], metric]
            /. metric -> CachedGCTensor[holder, metric]
            // ExpandPDToGCTensor[chart]
            // ContractGCTensors[holder]
        ,
            {a, -b, -c}
        ] // Fold[ChangeCovD[#1, PD, #2] &, #, subCDs] &
        // GCTensorHolderAction[holder, "PostCurvatureTensorCalculation"[chris], #] &;

        (* Riemann tensor *)
        holder /: HeldGCTensor[holder, riem] := ETensor[
            ChangeCurvature[riem[-a, -b, -c, -d], cd2, PD]
            /. {metric -> CachedGCTensor[holder, metric], chris -> CachedGCTensor[holder, chris]}
            // ExpandPDToGCTensor[chart]
            // ContractGCTensors[holder]
        ,
            {-a, -b, -c, -d}
        ] // RecoverSubRiemannTensors
        // ToCanonical[#, UseMetricOnVBundle -> None] &
        // Fold[ChristoffelToGradMetric, #, subMetrics] &
        // Fold[ChangeCovD[#1, PD, #2] &, #, subCDs] &
        // ContractMetric[#, subMetrics] &
        // GCTensorHolderAction[holder, "PostCurvatureTensorCalculation"[riem], #] &;

        (* Ricci tensor *)
        holder /: HeldGCTensor[holder, ricci] := ETensor[
            CachedGCTensor[holder, riem][-a, -c, -b, -d] CachedGCTensor[holder, metric][c, d] // ContractGCTensors[holder]
        ,
            {-a, -b}
        ] // ContractMetric[#, subMetrics] &
        // GCTensorHolderAction[holder, "PostCurvatureTensorCalculation"[ricci], #] &;

        (* Ricci scalar *)
        holder /: HeldGCTensor[holder, ricciScalar] := GCTensor[
            CachedGCTensor[holder, ricci][-a, -b] CachedGCTensor[holder, metric][a, b]
            // ContractGCTensors[holder]
            // Simplify
            // GCTensorHolderAction[holder, "PostCurvatureTensorCalculation"[ricciScalar], #] &
        ,
            {}
        ];
    ];
];
SyntaxInformation[AddCurvatureTensorsToHolder] = {"ArgumentsPattern" -> {_, _, _}};
Protect[AddCurvatureTensorsToHolder];

GetAllHeldTensors[holder_] := Union[
    Replace[#[[1]], {
        HoldPattern[Verbatim[HoldPattern][CachedGCTensor[holder, t_, _]]] :> t,
        HoldPattern[Verbatim[HoldPattern][HeldGCTensor[holder, t_]]] :> t,
        _ -> Nothing
    }] & /@ UpValues[holder]
];
SyntaxInformation[GetAllHeldTensors] = {"ArgumentsPattern" -> {_}};
Protect[GetAllHeldTensors];

GetAllHeldTensorRules[holder_] := # -> CachedGCTensor[holder, #] & /@ GetAllHeldTensors[holder];
SyntaxInformation[GetAllHeldTensorRules] = {"ArgumentsPattern" -> {_}};
Protect[GetAllHeldTensorRules];

ClearGCTensorHolderCache[holder_] := With[{
    t = Replace[#[[1]], {
        HoldPattern[Verbatim[HoldPattern][CachedGCTensor[holder, t_, inds_]]] /; !PersistentTensorCacheQ[holder, t] :> {t, inds},
        _ -> None
    }] & /@ UpValues[holder]
}, (holder /: CachedGCTensor[holder, #1, #2] =.) & @@@ t;];
SyntaxInformation[ClearGCTensorHolderCache] = {"ArgumentsPattern" -> {_}};
Protect[ClearGCTensorHolderCache];

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
    If[cperms != Range@Length@basis,
        MapIndexed[Function[{elem, inds},
            ETensorTranspose[elem, InversePermutation@Ordering[
                cperms[[
                    #[[1]] & /@ Position[inds - subMStartI, _?(# > 0 &), {1}]
                ]]
            ]]
        ], array, {Length@basis}] // Transpose[#, perms] &
    ,
        array
    ]
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
    MapThread[
        If[!ContractableChartsQ[#1, #2], Throw@Message[GCTensorContractTwo::icpbs, #1, #2]] &,
        {basis1[[#]] & /@ n1, basis2[[#]] & /@ n2}
    ];
    (* perm1 = With[{range = Range@d1}, InversePermutation@Join[Delete[range, Transpose@{n1}], range[[n1]]]]; *)
    perm1 = InversePermutation@Join[Delete[Range@d1, Transpose@{n1}], n1];
    perm2 = InversePermutation@Join[n2, Delete[Range@d2, Transpose@{n2}]];
    (* perm2 = With[{range = Range@d2}, InversePermutation@Join[range[[n2]], Delete[range, Transpose@{n2}]]]; *)
    If[$xDecompVerbose, Print["[GCTensorContractTwo] perm1 = ", perm1, ", perm2 = ", perm2]];
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
    perm = With[{range = Range@d}, InversePermutation@Join[Delete[range, {{n1}, {n2}}], range[[{n1, n2}]]]];
    If[$xDecompVerbose, Print["[GCTensorContrac] perm1 = ", perm]];
    GCTensorContractLast2[GCTensorTranspose[GCTensor[arr, basis], perm]]
];
SyntaxInformation[GCTensorContract] = {"ArgumentsPattern" -> {_, _, _}};
Protect[GCTensorContract];

GCTensorProduct[GCTensor[arr1_, basis1_], GCTensor[arr2_, basis2_]] := GCTensor[
    Map[Function[elem,
        Map[ETensorProduct[elem, #] &, arr2, {Length@basis2}]
    ], arr1, {Length@basis1}],
    Join[basis1, basis2]
];
SyntaxInformation[GCTensorProduct] = {"ArgumentsPattern" -> {_, _}};
Protect[GCTensorProduct];

GCTensorFixedContract[GCTensor[arr_, basis_], metric_GCTensor, n_] := GCTensorTranspose[
    GCTensorContractTwo[GCTensor[arr, basis], metric, {n}, {1}],
    MoveTo[Length@basis, n]
];
SyntaxInformation[GCTensorFixedContract] = {"ArgumentsPattern" -> {_, _, _}};
Protect[GCTensorFixedContract];

BasisETensorsFromSignedBasis[basis_Symbol] := Join[
    ToETensor@#[[1]] & /@ CoordsOfGChart@basis,
    ToETensor[#, {-1, 1}] & /@ SubManifoldsOfGChart@basis
];
BasisETensorsFromSignedBasis[-basis_Symbol] := Join[
    ToETensor@#[[2]] & /@ CoordsOfGChart@basis,
    ToETensor[#, {1, -1}] & /@ SubManifoldsOfGChart@basis
];
MoveIndicesOfVBundleToLeft[ETensor[expr_, inds_], vb_] := With[{
    pos = # & @@@ Position[vb === (VBundleOfIndex@#) & /@ inds, True, {1}]
}, ETensor[
    expr,
    Permute[inds, InversePermutation@Join[pos, Delete[Range@Length@inds, Transpose@{pos}]]]
]];
GCTensorToBasis[GCTensor[arr_, basis_]] := With[{
    bTensors = BasisETensorsFromSignedBasis /@ basis,
    clens = Length@CoordsOfGChart@UpBasis@# & /@ basis,
    vb = VBundleOfGChart@UpBasis@First@basis
}, MapIndexed[Function[{elem, indices},
    ETensorContractTwo[
        MoveIndicesOfVBundleToLeft[ETensorProduct @@ MapThread[Part, {bTensors, indices}], vb],
        elem,
        Length@Position[Thread[indices > clens], True, {1}]
    ]
], arr, {Length@basis}] // Total[#, Length@basis] &];
SyntaxInformation[GCTensorToBasis] = {"ArgumentsPattern" -> {_}};
Protect[GCTensorToBasis];

AIndexToPattern[a_Symbol] := Pattern[a, Blank[]];
AIndexToPattern[-a_Symbol] := -Pattern[a, Blank[]];
CoordParamsOfGChart[chart_] := #3 & @@@ CoordsOfGChart@UpBasis@chart;
GCTensorElementToRule[elem_, indices_, cparams_, head_] := (
    (HoldPattern[head[##]] & @@ MapThread[Part, {cparams, indices}]) -> elem
) /; And @@ Thread[indices <= Length /@ cparams];
GCTensorElementToRule[ETensor[expr_, inds_], indices_, cparams_, head_] := With[{
    indList = First@Fold[With[{
        last = #1[[1]],
        cursor = #1[[2]],
        i = #2[[1]],
        cparam = #2[[2]]
    },
        If[i > Length@cparam,
            {Append[last, AIndexToPattern@inds[[cursor]]], cursor + 1},
            {Append[last, cparam[[i]]], cursor}
        ]
    ] &, {{}, 1}, Thread@{indices, cparams}],
    dummies = List @@ Union[UpIndex /@ FindDummyIndices@expr]
},
    (HoldPattern[head[##]] & @@ indList) :> Module[dummies, expr]
];
SortComponentRules[_ -> 0] = 0;
SortComponentRules[_ :> Module[{}, 0]] = 1;
SortComponentRules[_] = 2;
SetAttributes[SortComponentRules, HoldFirst];

Options[GCTensorToComponentRules] = {
    HeadOfCoordinateIndex -> LI,
    ExcludeZeros -> False
};
GCTensorToComponentRules[GCTensor[arr_, basis_], head_, opt: OptionsPattern[]] := Module[
    {li, cparams, dim, res},
    li = OptionValue[HeadOfCoordinateIndex];
    cparams = MapThread[If[MatchQ[#1, _Symbol], li /@ #2, -(li /@ #2)] &, {basis, CoordParamsOfGChart /@ basis}];
    dim = Length@basis;
    res = MapIndexed[Function[{elem, indices},
        GCTensorElementToRule[elem, indices, cparams, head]
    ], arr, {dim}] // If[dim > 1, Flatten[#, dim - 1], #] &;
    If[OptionValue[ExcludeZeros],
        DeleteCases[res, Alternatives[_ -> 0, _ :> Module[{}, 0]]]
    ,
        SortBy[res, SortComponentRules]
    ]
];
SyntaxInformation[GCTensorToComponentRules] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
Protect[GCTensorToComponentRules];

DefGCTensorMapFunc[funcs__] := Function[func,
    GCTensor /: func[GCTensor[arr_, basis_], args___] := GCTensor[Map[func[#, args] &, arr, {Length@basis}], basis];
    GCTensor /: func[GCTensor[arr_, basis_][inds__], args___] := GCTensor[Map[func[#, args] &, arr, {Length@basis}], basis][inds];
] /@ {funcs};

DefGCTensorMapFunc[
    ScreenDollarIndices,
    ToCanonical,
    ContractMetric,
    ChangeCovD,
    ChangeCurvature,
    ChristoffelToRiemann,
    ChristoffelToGradMetric
];
GCTensor[e_, {}][] := Scalar[e];
GCTensor /: ParamD[args__]@GCTensor[arr_, basis] := GCTensor[Map[ParamD[args], arr, {Length@basis}], basis];
GCTensor /: ETensor[t_GCTensor[inds__], {inds2__}] := With[{
    perm = PermutationProduct[InversePermutation@Ordering@{inds}, Ordering@{inds2}]
}, GCTensorTranspose[t, perm]] /; Sort@{inds} === Sort@{inds2};

GCTensor /: GCTensor[arr1_, basis_] + GCTensor[arr2_, basis_] := GCTensor[arr1 + arr2, basis];
GCTensor /: x_?xTools`xTension`Private`IndexedScalarQ * GCTensor[arr_, basis_] := GCTensor[arr * x, basis];
GCTensor /: GCTensor[arr1_, basis1_][inds1__] + GCTensor[arr2_, basis2_][inds2__] := With[{
    perm = PermutationProduct[InversePermutation@Ordering@{inds2}, Ordering@{inds1}]
},
    GCTensor[arr1 + GCTensorDataTranspose[arr2, basis2, perm], basis1][inds1]
] /; Sort@{inds1} === Sort@{inds2} && Length@basis1 === Length@basis2;

GCTensor /: x_?xAct`xTensor`Private`NonIndexedScalarQ * GCTensor[arr_, basis_][inds__] := GCTensor[arr * x, basis][inds];

SyntaxInformation[GCTensor] = {"ArgumentsPattern" -> {_, _}};
Protect[GCTensor];

ValidateGCTensor[GCTensor[arr_, basis_]] := With[{
    subMs = SubManifoldsOfGChart@UpBasis@# & /@ basis,
    pms = If[MatchQ[#, _Symbol], 1, -1] & /@ basis
}];
SyntaxInformation[ValidateGCTensor] = {"ArgumentsPattern" -> {_}};
Protect[ValidateGCTensor];

SignedVBundleOfIndex[id_Symbol] := VBundleOfIndex@id;
SignedVBundleOfIndex[-id_Symbol] := -VBundleOfIndex@id;

(* CreateGCTensorOneComponent[{coords__} -> expr_, charts_List] := ; *)
PartOrNone[expr_, n_Integer /; n > 0] := expr[[n]];
PartOrNone[_, n_Integer /; n <= 0] = Nothing;
PartSpecFromComponentSpec[spec_List, params_, subVBundles_] := MapThread[
    {s, p, v} |-> With[{
        ppos = FirstPosition[p, s]
    }, If[!MatchQ[ppos, _Missing],
        ppos[[1]]
    ,
        Length@p + FirstPosition[subVBundles, SignedVBundleOfIndex@s][[1]]
    ]]
, {spec, params, subVBundles}];
InitGCTensorElement[clens_, subVBundles_][indices__] := ZeroETensor@MapThread[PartOrNone, {subVBundles, {indices} - clens}];
AddGCTensorElement[params_, subVBundles_][arr_, spec_List -> expr_] := Module[
    {inds, ainds, elem},
    inds = PartSpecFromComponentSpec[spec, params, subVBundles];
    clens = Length /@ params;
    ainds = MapThread[If[#2 > 0, #1, Nothing] &, {spec, inds - clens}];
    With[{
        ainds2 = List @@ FindFreeIndices@expr
    }, If[ainds =!= ainds2, Throw@Message[CreateGCTensor::unmatched, ainds, ainds2]]];
    elem = If[Length@ainds > 0, ETensor[expr, ainds], expr];
    ReplacePart[arr, inds -> elem]
];
CreateGCTensor::unmatched = "Indices `1` and `2` don't match.";
CreateGCTensor[{components__}, charts_] := Module[
    {dim, arrDim, params, subVBundles, ret},
    dim = Length@charts;
    params = MapThread[If[MatchQ[#1, _Symbol], #2, -#2] &, {charts, CoordParamsOfGChart /@ charts}];
    subVBundles = SubvbundlesOfGChart /@ charts;
    ret = Fold[
        AddGCTensorElement[params, subVBundles],
        Array[InitGCTensorElement[Length /@ params, subVBundles], (Length /@ params) + (Length /@ subVBundles)],
        {components}
    ];
    GCTensor[ret, charts]
];
SyntaxInformation[CreateGCTensor] = {"ArgumentsPattern" -> {_, _}};
Protect[CreateGCTensor];

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

ProductTwoIndexedGCTensors[t1_GCTensor[inds1__], t2_GCTensor[inds2]] := GCTensorProduct[t1, t2] @@ Join[{inds1}, {inds2}];

ExecuteContraction[{_, _, idx_, idx_}, l_List] := {
    ContractOneIndexedGCTensors[l[[idx]]],
    Delete[l, idx]
};
ExecuteContraction[{_, c_, idx1_, idx2_}, l_List] := With[{
    t1 = l[[idx1]],
    t2 = l[[idx2]]
}, {
    If[c > 0, ContractTwoIndexedGCTensors[t1, t2], ProductTwoIndexedGCTensors[t1, t2]],
    Delete[l, {{idx1}, {idx2}}]
}] /; idx1 != idx2;
OptimizedGCTensorContraction[{}] = 1;
OptimizedGCTensorContraction[l_List] := Module[
    {tensors, contractions, contractionCounts, selected, res, l2},
    tensors = MapIndexed[{List @@ #1, #2[[1]]} &, Cases[l, GCTensor[__][__]]];
    If[Length@tensors > 0,
        contractions = Flatten[
            Outer[{
                Length@xAct`xTensor`Private`DropPairs[#1[[1]], #2[[1]]],
                Length@xAct`xTensor`Private`TakePairs[#1[[1]], #2[[1]]],
                #1[[2]],
                #2[[2]]
            } &, tensors, tensors, 1],
            1
        ];
        contractionCounts = #[[1]] & /@ contractions;
        selected = contractions[[FirstPosition[contractionCounts, Min@contractionCounts][[1]]]];
        If[selected[[2]] > 0 || selected[[3]] != selected[[4]],
            If[$xDecompVerbose, Print["[OptimizedGCTensorContraction]", "selected: ", tensors[[selected[[3]], 1]], tensors[[selected[[4]], 1]]]];
            {res, l2} = ExecuteContraction[selected, l];
            OptimizedGCTensorContraction@Append[l2, res]
        ,
            Times @@ l
        ]
    ,
        Times @@ l
    ]
];

ContractGCTensors[covd_][expr_] := ContractGCTensors[expr, covd];
ContractGCTensors[expr_Plus, covd_] := ContractGCTensors[covd] /@ expr;
ContractNonTensors[e: GCTensor[__][__], _] := e;
ContractNonTensors[expr_, covd_] := ContractGCTensors[expr, covd];
ContractGCTensors[expr_Times, covd_] := With[
    {l = ContractNonTensors[#, covd] & /@ List @@ expr},
    (Times @@ DeleteCases[l, GCTensor[__][__]]) * OptimizedGCTensorContraction[Cases[l, GCTensor[__][__]]]
];
ContractGCTensors[expr: GCTensor[__][__], _] := OptimizedGCTensorContraction@{expr};
ContractGCTensors[expr_, _] := expr;
SyntaxInformation[ContractGCTensors] = {"ArgumentsPattern" -> {_, _.}};
Protect[ContractGCTensors];

RecoverOneRiemannETensor[ETensor[expr_, {-a_Symbol, -b_Symbol, -c_Symbol, -d_Symbol}]] := Module[
    {vb, cd, chris, riem, rep, rule},
    vb = VBundleOfIndex@a;
    cd = CovDOfMetric@First@MetricsOfVBundle@vb;
    chris = Christoffel[cd][a, -b, -c][[0]];
    riem = Riemann@cd;
    rep = ChangeCurvature[riem[-a, -b, -c, d], cd] - riem[-a, -b, -c, d] + PD[-a][chris[d, -b, -c]];
    rule = {
        PD[-a]@chris[d_Symbol, -b, -c] :> Module[#1, #2],
        PD[-a]@chris[d_Symbol, -c, -b] :> Module[#1, #2]
    } &[Union[List @@ (UpIndex /@ FindDummyIndices@Evaluate@rep)], rep];
    ETensor[expr /. rule, {-a, -b, -c, -d}]
];
RecoverSubRiemannTensors[GCTensor[arr_, {-chart_?GChartQ, -chart_, -chart_, -chart_}]] := Module[
    {clen, mlen},
    clen = Length@CoordsOfGChart@chart;
    mlen = Length@SubManifoldsOfGChart@chart;
    MapAt[RecoverOneRiemannETensor, arr, {#, #, #, #} & /@ Range[clen + 1, clen + mlen]] // GCTensor[#, -{chart, chart, chart, chart}] &
];
SyntaxInformation[RecoverSubRiemannTensors] = {"ArgumentsPattern" -> {_}};
Protect[RecoverSubRiemannTensors];

On[RuleDelayed::rhs];

End[];

EndPackage[];