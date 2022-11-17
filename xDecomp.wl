BeginPackage["xTools`xDecomp`", {"xAct`xCore`", "xAct`xTensor`", "xTools`xTension`"}];

(Unprotect[#]; ClearAll[#];) & /@ Names@{$Context <> "*", $Context <> "Private`*"};

GChartQ::usage = "GChartQ[chart] gives True if chart is a generalized chart.";
CoordsOfGChart::usage = "CoordsOfGChart[chart] gives single-variable coordinate list of the generalized chart, each element is in the form {e, ed, param}."
SubManifoldsOfGChart::usage = "SubManifoldsOfGChart[chart] gives a list of basis tensors of the generalized chart.";
MetricOfGChart::usage = "MetricOfGChart[chart] gives the metric of generalized chart, if defined."
GChartMetricDecompose::usage = "GChartMetricDecompose[chart, metric[-a, -b]] decomposes the metric using the defined metric decomposition rule.";
VBundleOfGChart::usage = "VBundleOfGChart[chart] gives the base VBundle of the generalized chart.";
SubVBundlesOfGChart::usage = "SubVBundlesOfGChart[chart] gives the subvbundles of sub manifolds of the generalized chart.";
IndexDimensionOfGChart::usage = "IndexDimensionOfGChart[chart] gives the index dimension, i.e., the dimension where sub manifolds are treated as dim 1, of the chart.";
GetBasisETensorsOfGChart::usage = "GetBasisETensorsOfGChart[chart] gives a list of basis of chart in ETensors form.";

DefGBasis::usage = "DefGBasis[chart, {{e, ed, param}...}, {eT...}] defines a generalized chart with coodinate basis and sub manifolds."
UndefGBasis::usage = "UndefGBasis[chart] undefines a generalized chart.";
SetGChartMetric::usage = "SetGChartMetric[chart, metric, metricValue, metricInvValue] defines metric decomposition rule for metric.";

GetAllHeldTensors::usage = "GetAllCachedTensors[holder] gives a list of tensors that has a definition with the holder.";
GetAllHeldTensorRules::usage = "GetAllHeldTensorRules[holder] gives a rule list to replace all tensors defined in holder as the values.";
ClearGCTensorHolderCache::usage = "ClearGCTensorHolderCache[holder] clears the tensor cache.";

GCTensorHolderQ::usage = "GCTensorHolderQ[covd] gives True for GCovD.";

SetHeldMetric::usage = "SetHeldMetric[holder, metric, metric, invMetric] sets the metric of the holder.";
AddCurvatureTensorsToHolder::usage = "AddCurvatureTensorsToHolder[holder, chart, christoffelTensor] defines curvature tensors to the holder, using the given Christoffel tensor.";
SortGChartParamD::usage = "SortGChartParamD[expr, chart] or SortGChartParamD[chart][expr] makes all coordinate ParamDs of chart in expr insider PDs.";

HeldGCTensor::usage = "HeldGCTensor[holder, tensor] gives the tensor provider of the tensor.";
CachedGCTensor::usage = "CachedGCTensor[holder, tensor, inds] gives the cached tensor with specific risen/lowered indices.";

GCTensorHolderAction::usage = "GCTensorHolderAction[holder, tag, expr] defines some extra action applies in various computation stages. Typically used to simplify expressions.";
GCTensorHolderDefaultAction::usage = "GCTensorHolderDefaultAction[holder, tag, expr] is the default value for GCTensorHolderAction[holder, tag, expr].";
GCTensorHolderDAUseMetricVB::usage = "GCTensorHolderDAUseMetricVB[holder] is the UseMetricOnVBundle option passed to ToCanonical by GCTensorHolderDefaultAction[holder]. Default is All.";
PostCurvatureTensorCalculation::usage = "PostCurvatureTensorCalculation[tensor] is an action tag for GCTensorHolderAction that's executed after the calculation of every curvature related tensors.";
PostChangeIndex::usage = "PostChangeIndex[tensor, indices, changedInd] is an action tag for GCTensorHolderAction that's executed after changing index for a cached tensor.";

PDGChart::usage = "PDGChart[chart, -a] represents the PD operator of the GChart chart.";
ExpandPDToBasis::usage = "ExpandPDToBasis[expr, chart] or ExpandPDToBasis[chart][expr] expands all partial derivatives in expr in terms of basis the coordinates and sub manifolds of chart.";
ExpandMetric::usage = "ExpandMetric[expr, chart] or ExpandMetric[chart][expr] replaces all metric in expr with the decomposed one defined in chart."
ExpandPDToGCTensor::usage = "ExpandPDToGCTensor[expr, chart] or ExpandPDToGCTensor[chart][expr] replaces all PDs to GCTensors in expr.";

GCTensor::usage = "GCTensor[values, charts] represents a generalized CTensor.";
GCTensorChangeBasis::usage = "GCTensorChangeBasis[GCTensor[...], n, basis] changes the basis of n-th axis to the specified basis.";
EnsureGCTensorBasis::usage = "EnsureGCTensorBasis[GCTensor[...], basis] changes the basis of the GCTensor of any of the axes doesn't match the given one.";
GCTensorTranspose::usage = "GCTensorTranspose[GCTensor[...], perms] performs the tensor transpose of GCTensor.";
GCTensorContract::usage = "GCTensorContract[T, n1, n2] contracts the n1 axis with n2 axis of T.";
GCTensorContractTwo::usage = "GCTensorContractTwo[T1, T2, n1, n2] computes the tensor contraction between two tensors, in which n1-th axis of T1 and n2-th axis of T2.";
GCTensorProduct::usage = "GCTensorProduct[T1, T2] computes the tensor product of T1 and T2.";
GCTensorFixedContract::usage = "GCTensorFixedContract[T, T2, n] contracts n-th axis of T with the symmetric tensor T2, with the resulting axis stay at the original position. Typically used to change index.";
GCTensorToBasis::usage = "GCTensorToBasis[expr] converts GCTensor expressions to basis.";
ZeroGCTensorQ::usage = "ZeroGCTensorQ[T] gives True of GCTensor T is zero.";
PostETensorContract::usage = "PostETensorContract is an option of GCTensorContractTwo, GCTensorContract, ContractGCTensors that's called after each ETensor contraction.";

GCTensorToComponentRules::usage = "GCTensorToComponentRules[T, head] converts GCTensor to a list of component rules.";
HeadOfCoordinateIndex::usage = "HeadOfCoordinateIndex is an option of GCTensorToComponentRules that defines the head used for coordinate indices. Default is LI.";
ExcludeZeros::usage = "ExcludeZeros is a boolean option of GCTensorToComponentRules, default to True. When True, only non-zero components are listed.";

ContractGCTensors::usage = "ContractGCTensors[expr, covd] performs contractions of GCTensor's, use metric if needed.";
ValidateGCTensor::usage = "ValidateGCTensor[GCTensor[...]] fixes some ";
ZeroGCTensor::usage = "ZeroGCTensor[charts] creats a GCTensor with zero components.";
CreateGCTensor::usage = "CreateGCTensor[{{r1, r2, ...} -> expr}, charts] is a convenient function for creating GCTensors.";
GCTensorPDGrad::usage = "GCTensorPDGrad[expr, chart] calculates the partial derivative of the tensor expression in chart.";
GCTensorPDDiv::usage = "GCTensorPDDiv[expr, n, chart] acts the PD operator on expr and contracts it with n-th axis.";

RecoverSubRiemannTensors::usage = "RecoverSubRiemannTensors[GCTensor[..., {-c1, -c2, -c3, -c4}]] replaces Christoffel tensors of sub manifolds with Riemann tensors.";

$xDecompVerbose = False;
$xDecompVerbose::usage = "$xDecompVerbose is a global boolean variable, when set to True, debug messages wil be printed.";

$xDecompDebugPrint = Print;
$xDecompDebugPrint::usage = "$xDecompDebugPrint is a global hook variable for printing debug messages.";

Begin["`Private`"];

Off[RuleDelayed::rhs];

DebugPrint[msg___] := If[$xDecompVerbose, $xDecompDebugPrint["[xDecomp] ", msg]];
SetAttributes[DebugPrint, HoldAll];

GChartQ[_] = False;
CoordsOfGChart[_] = None;
CoordsOfGChart[-a_Symbol] := CoordsOfGChart[a];
SubManifoldsOfGChart[_] = None;
MetricOfGChart[_] = None;
GChartParamQ[_] = False;
GChartMetricDecompose[_, expr_] := expr;
VBundleOfGChart[chart_] := First@SlotsOfTensor@CoordsOfGChart[chart][[1, 1]];
SubVBundlesOfGChart[chart_Symbol] := SlotsOfTensor[#][[1]] & /@ SubManifoldsOfGChart[chart];
SubVBundlesOfGChart[-chart_Symbol] := -SubVBundlesOfGChart[chart];
SyntaxInformation[GChartQ] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[CoordsOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SubManifoldsOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[MetricOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[GChartMetricDecompose] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[VBundleOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SubVBundlesOfGChart] = {"ArgumentsPattern" -> {_}};

IndexDimensionOfGChart[chart_?GChartQ] := Length[CoordsOfGChart[chart]] + Length[SubManifoldsOfGChart[chart]];
SyntaxInformation[IndexDimensionOfGChart] = {"ArgumentsPattern" -> {_}};

ChangeGChartSign[a_?GChartQ] := -a;
ChangeGChartSign[-a_?GChartQ] := a;
SignOfGChart[a_?GChartQ] = 1;
SignOfGChart[-a_?GChartQ] = -1;

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

UndefGBasis[l_List] := UndefGBasis /@ l;
UndefGBasis[chart_Symbol?GChartQ] := (
    Apply[{e, ed, param} |-> With[{
        Q = Symbol[ToString[First@SlotsOfTensor@e] <> "`Q"]
    },
        param /: PD[-a_Symbol?Q][param] =.;
        UndefTensor /@ {e, ed};
    ], CoordsOfGChart@chart, {1}];
    UndefTensor /@ SubManifoldsOfGChart@chart;
    Remove@chart;
);
SyntaxInformation[UndefGBasis] = {"ArgumentsPattern" -> {___}};

IndexToPattern[idx_Symbol] := Pattern[idx, Blank[]];
IndexToPattern[-idx_Symbol] := -Pattern[idx, Blank[]];
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

GetBasisETensorsOfGChart[chart_?GChartQ] := With[{
    ai = First@GetIndicesOfVBundle[VBundleOfGChart@chart, 1]
}, Join[
    ETensor[#[[1]][ai], {ai}] & /@ CoordsOfGChart@chart,
    With[{bi = First@GetIndicesOfVBundle[First@SlotsOfTensor@#, 1]},
        ETensor[#[-bi, ai]]
    ] & /@ SubManifoldsOfGChart@chart
]];
GetBasisETensorsOfGChart[-chart_?GChartQ] := With[{
    ai = First@GetIndicesOfVBundle[VBundleOfGChart@chart, 1]
}, Join[
    ETensor[#[[2]][-ai], {-ai}] & /@ CoordsOfGChart@chart,
    With[{bi = First@GetIndicesOfVBundle[First@SlotsOfTensor@#, 1]},
        ETensor[#[bi, -ai]]
    ] & /@ SubManifoldsOfGChart@chart
]];
SyntaxInformation[GetBasisETensorsOfGChart] = {"ArgumentsPattern" -> {_}};

GCTensorHolderQ[_] = None;
SyntaxInformation[GCTensorHolderQ] = {"ArgumentsPattern" -> {_}};

UpDownIndexToNumber[_Symbol] = 1;
UpDownIndexToNumber[-_Symbol] = -1;
UpDownIndexToNumber[_] = 0;

HeldGCTensor::undef = "tensor provider for `1` of holder `2` is not defined.";
HeldGCTensor[holder_, tensor_] := Throw@Message[HeldGCTensor::undef, tensor, holder];
SyntaxInformation[HeldGCTensor] = {"ArgumentsPattern" -> {_, _}};

GCTensorHolderDAUseMetricVB[_] = All;
SyntaxInformation[GCTensorHolderDAUseMetricVB] = {"ArgumentsPattern" -> {_}};

GCTensorHolderAction[holder_, tag_, expr_] := GCTensorHolderDefaultAction[holder, tag, expr];
GCTensorHolderDefaultAction[holder, _, e_] := e;
GCTensorHolderDefaultAction[holder_, PostChangeIndex[tensor_, inds_, i_], expr_] := Simplify@ToCanonical[ContractMetric[expr], UseMetricOnVBundle -> GCTensorHolderDAUseMetricVB@holder];
GCTensorHolderDefaultAction[holder_, PostChangeIndex[tensor_, inds_, i_], expr: GCTensor[arr_, basis_]] := Simplify@ToCanonical[ContractMetric[expr, First@MetricsOfVBundle@# & /@ SubVBundlesOfGChart@UpBasis@basis[[i]]], UseMetricOnVBundle -> GCTensorHolderDAUseMetricVB@holder];
GCTensorHolderDefaultAction[holder_, PostCurvatureTensorCalculation[tensor_], expr_] := Simplify@ToCanonical[expr, UseMetricOnVBundle -> None];
SyntaxInformation[GCTensorHolderAction] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[GCTensorHolderDefaultAction] = {"ArgumentsPattern" -> {_, _, _}};

SyntaxInformation[PostChangeIndex] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[PostCurvatureTensorCalculation] = {"ArgumentsPattern" -> {_}};

CachedGCTensor[holder_Symbol, tensor_][inds___] := CachedGCTensor[holder, tensor, UpDownIndexToNumber /@ {inds}][inds];
CachedGCTensor[holder_Symbol, tensor_, inds_] := (
    DebugPrint["[CachedGCTensor] calculating tensor ", tensor];
    holder /: CachedGCTensor[holder, tensor, inds] = HeldGCTensor[holder, tensor]
) /; inds === UpDownIndexToNumber /@ SlotsOfTensor@tensor;
CachedGCTensor[holder_Symbol, tensor_, inds_] := With[{
    firstUnmatched = First@FirstPosition[MapThread[#1 =!= UpDownIndexToNumber@#2 &, {inds, SlotsOfTensor@tensor}], True],
    slots = SlotsOfTensor@tensor
},
    holder /: CachedGCTensor[holder, tensor, inds] = With[{
        tensor2 = CachedGCTensor[holder, tensor, MapAt[-# &, inds, firstUnmatched]]
    },
        DebugPrint["[CachedGCTensor] ", "calculating index change ", tensor, inds, " from ", MapAt[-# &, inds, firstUnmatched]];
        GCTensorFixedContract[
            tensor2,
            CachedGCTensor[holder, First@MetricsOfVBundle@slots[[firstUnmatched]], {#, #} &[inds[[firstUnmatched]]]],
            firstUnmatched
        ] // GCTensorHolderAction[holder, PostChangeIndex[tensor, inds, firstUnmatched], #] &
    ]
] /; inds =!= UpDownIndexToNumber /@ SlotsOfTensor@tensor;
SyntaxInformation[CachedGCTensor] = {"ArgumentsPattern" -> {_, _, _.}};

PersistentTensorCacheQ[_, _] = False;
SetHeldMetric[holder_Symbol, metric_, metrici_GCTensor, invMetric_GCTensor] := (
    holder /: CachedGCTensor[holder, metric, {-1, -1}] = metrici;
    holder /: CachedGCTensor[holder, metric, {1, 1}] = invMetric;
    holder /: PersistentTensorCacheQ[holder, metric] = True;
);
SyntaxInformation[SetHeldMetric] = {"ArgumentsPattern" -> {_, _, _, _}};

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
        subCDs = CovDOfMetric@First@MetricsOfVBundle@# & /@ SubVBundlesOfGChart@chart,
        subMetrics = First@MetricsOfVBundle@# & /@ SubVBundlesOfGChart@chart,
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
        // GCTensorHolderAction[holder, PostCurvatureTensorCalculation[chris], #] &;

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
        // SortGChartParamD[chart]
        // Fold[ChangeCovD[#1, PD, #2] &, #, subCDs] &
        // ContractMetric[#, subMetrics] &
        // GCTensorHolderAction[holder, PostCurvatureTensorCalculation[riem], #] &;

        (* Ricci tensor *)
        holder /: HeldGCTensor[holder, ricci] := ETensor[
            CachedGCTensor[holder, riem][-a, -c, -b, -d] CachedGCTensor[holder, metric][c, d] // ContractGCTensors[holder]
        ,
            {-a, -b}
        ] // ContractMetric[#, subMetrics] &
        // GCTensorHolderAction[holder, PostCurvatureTensorCalculation[ricci], #] &;

        (* Ricci scalar *)
        holder /: HeldGCTensor[holder, ricciScalar] := GCTensor[
            CachedGCTensor[holder, ricci][-a, -b] CachedGCTensor[holder, metric][a, b]
            // ContractGCTensors[holder]
            // Simplify
            // GCTensorHolderAction[holder, PostCurvatureTensorCalculation[ricciScalar], #] &
        ,
            {}
        ];
    ];
];
SyntaxInformation[AddCurvatureTensorsToHolder] = {"ArgumentsPattern" -> {_, _, _}};

SortGChartParamD[chart_][expr_] := SortGChartParamD[expr, chart];
SortGChartParamD[expr_, chart_] := With[{
    rules = (ParamD[l___, #3, r___]@PD[-a_Symbol]@A_ :> ParamD[l, r]@PD[-a]@ParamD[#3]@A) & @@@ CoordsOfGChart@chart
}, expr //. rules];
SyntaxInformation[SortGChartParamD] = {"ArgumentsPattern" -> {_, _.}};

GetAllHeldTensors[holder_] := Union[
    Replace[#[[1]], {
        HoldPattern[Verbatim[HoldPattern][CachedGCTensor[holder, t_, _]]] :> t,
        HoldPattern[Verbatim[HoldPattern][HeldGCTensor[holder, t_]]] :> t,
        _ -> Nothing
    }] & /@ UpValues[holder]
];
SyntaxInformation[GetAllHeldTensors] = {"ArgumentsPattern" -> {_}};

GetAllHeldTensorRules[holder_] := # -> CachedGCTensor[holder, #] & /@ GetAllHeldTensors[holder];
SyntaxInformation[GetAllHeldTensorRules] = {"ArgumentsPattern" -> {_}};

ClearGCTensorHolderCache[holder_] := With[{
    t = Replace[#[[1]], {
        HoldPattern[Verbatim[HoldPattern][CachedGCTensor[holder, t_, inds_]]] /; !PersistentTensorCacheQ[holder, t] :> {t, inds},
        _ -> None
    }] & /@ UpValues[holder]
}, (holder /: CachedGCTensor[holder, #1, #2] =.) & @@@ t;];
SyntaxInformation[ClearGCTensorHolderCache] = {"ArgumentsPattern" -> {_}};

PDGChart[chart_?GChartQ, -a_Symbol][t_GCTensor[inds__]] := GCTensorPDGrad[t, chart][inds, -a];
SyntaxInformation[PDGChart] = {"ArgumentsPattern" -> {_, _}};

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

ExpandPDToGCTensor[chart_][expr_] := ExpandPDToGCTensor[expr, chart];
ExpandPDToGCTensor[expr_, chart_?GChartQ] := With[{
    Q = Symbol[ToString[VBundleOfGChart@chart] <> "`Q"]
}, expr /. PD[-a_Symbol?Q] -> PDGChart[chart, -a]];
SyntaxInformation[ExpandPDToGCTensor] = {"ArgumentsPattern" -> {_, _.}};

ExpandMetric[chart_?GChartQ][expr_] := ExpandMetric[expr, chart];
ExpandMetric[expr_, chart_?GChartQ] := With[{metric = MetricOfGChart@chart},
   expr /. metric[inds__] :> GChartMetricDecompose[chart, metric[inds]]
];
SyntaxInformation[ExpandMetric] = {"ArgumentsPattern" -> {_, _.}};

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

GCTensorProduct[GCTensor[arr1_, basis1_], GCTensor[arr2_, basis2_]] := GCTensor[
    Map[Function[elem,
        Map[ETensorProduct[elem, #] &, arr2, {Length@basis2}]
    ], arr1, {Length@basis1}],
    Join[basis1, basis2]
];
SyntaxInformation[GCTensorProduct] = {"ArgumentsPattern" -> {_, _}};

GetBasisChangeMatrix[basisTo_, basisFrom_] := Outer[
    ETensorContractTwo[#1, #2, {-1}, {-1}] &,
    GetBasisETensorsOfGChart[basisTo],
    GetBasisETensorsOfGChart[ChangeGChartSign@basisFrom]
];

Null@ContractPair;

GCTensorDotInner[] := GCTensorDotInner[Identity];
GCTensorDotInner[action_][lhs_, ContractPair[rhs_, len_]] := With[{r1 = Range[-len, -1], r2 = Range[len]},
    ETensorContractTwo[lhs, rhs, r1, r2] // If[len > 0, action@#, #] &
];

GCTensorChangeBasis::invalid = "Cannot change basis with different signs `1` to `2`.";
GCTensorChangeBasis[GCTensor[arr_, basis_], n_Integer, basis2_] := With[{
    mat = GetBasisChangeMatrix[basis2, basis[[n]]],
    clen = Length@CoordsOfGChart@basis2
},
    If[SignOfGChart[basis2] != SignOfGChart@basis[[n]], Throw@Message[GCTensorChangeBasis::invalid, basis[[n]], basis2]];
    Map[Function[elem, Inner[
        GCTensorDotInner[],
        mat,
        MapIndexed[ContractPair[#1, If[#2[[1]] > clen, 1, 0]] &, elem],
        Plus
    ]], arr, {n - 1}] // GCTensor[#, ReplacePart[basis, n -> basis2]] &
];
SyntaxInformation[GCTensorChangeBasis] = {"ArgumentsPattern" -> {_, _}};

EnsureGCTensorBasis[GCTensor[arr_, basis_], basis2_List] := Fold[
    If[#2[[3]] =!= None && #2[[2]] =!= #2[[3]],
        GCTensorChangeBasis[#1, #2[[1]], #2[[3]]],
        #1
    ] &,
    GCTensor[arr, basis],
    Thread[{Range@Length@basis, basis, basis2}]
];
SyntaxInformation[EnsureGCTensorBasis] = {"ArgumentsPattern" -> {_, _}};

UpDownChartsQ[c_?GChartQ, -c_?GChartQ] = True;
UpDownChartsQ[-c_?GChartQ, c_?GChartQ] = True;
UpDownChartsQ[_, _] = False;

GCTensorDotSameBasis[GCTensor[arr1_, basis1_], GCTensor[arr2_, basis2_], len_, action_] := Module[
    {coordLengths, a1, a2},
    coordLengths = Length@CoordsOfGChart@UpBasis@# & /@ basis1[[Length[basis1] - len + 1 ;;]];
    a1 = If[len > 1, Map[Flatten, arr1, {Length[basis1] - len}], arr1];
    a2 = MapIndexed[Function[{elem, indices},
        Map[ContractPair[
            #,
            Length@Cases[Thread[indices > coordLengths], True]
        ] &, elem, {Length@basis2 - len}]
    ], arr2, {len}] // If[len > 1, Flatten[#, len - 1], #] &;
    GCTensor[
        Inner[GCTensorDotInner[action], a1, a2, Plus],
        Join[
            basis1[[;; Length[basis1] - len]],
            basis2[[len + 1 ;;]]
        ]
    ]
];

GCTensorContractTwo::icpbs = "cannot contract incompatible basis `1` and `2`.";
GCTensorContractTwo::argrx = "Expected 4 or more arguments, `1` found.";
Options[GCTensorContractTwo] = {PostETensorContract -> Identity};
GCTensorContractTwo[GCTensor[arr1_, basis1_], GCTensor[arr2_, basis2_], n1_List, n2_List, opt: OptionsPattern[]] := Module[
    {d1, d2, perm1, perm2, trans1, trans2},
    d1 = Length@basis1;
    d2 = Length@basis2;
    MapThread[
        If[!UpDownChartsQ[#1, #2], Throw@Message[GCTensorContractTwo::icpbs, #1, #2]] &,
        {basis1[[#]] & /@ n1, basis2[[#]] & /@ n2}
    ];
    perm1 = InversePermutation@Join[Delete[Range@d1, Transpose@{n1}], n1];
    perm2 = InversePermutation@Join[n2, Delete[Range@d2, Transpose@{n2}]];
    DebugPrint["[GCTensorContractTwo] perm1 = ", perm1, ", perm2 = ", perm2];
    trans1 = GCTensorTranspose[GCTensor[arr1, basis1], perm1];
    trans2 = GCTensorTranspose[GCTensor[arr2, basis2], perm2];
    trans2 = EnsureGCTensorBasis[trans2, Join[ChangeGChartSign /@ basis1[[n1]], ConstantArray[None, Length@basis2 - Length@n1]]];
    GCTensorDotSameBasis[trans1, trans2, Length@n1, OptionValue[PostETensorContract]]
] /; Length@n1 > 0;
GCTensorContractTwo[t1_GCTensor, t2_GCTensor, {}, {}, opt: OptionsPattern[]] := GCTensorProduct[t1, t2];
GCTensorContractTwo[args___] := Throw@Message[GCTensorContractTwo::argrx, Length@{args}];
SyntaxInformation[GCTensorContractTwo] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

GCTensorContractLast2SameBasis[GCTensor[arr_, basis_], action_] := With[{
    ccount = Length@CoordsOfGChart@UpBasis@basis[[-1]]
},
    Map[Function[elem, With[{diag = Diagonal@elem},
        Total[diag[[;; ccount]]] + Total[action@ETensorContract[#, -2, -1] & /@ diag[[ccount + 1 ;;]]]
    ]], arr, {Length@basis - 2}] // GCTensor[#, basis[[;; Length@basis - 2]]] &
];

GCTensorContract::icpbs = "cannot contract incompatible basis `1` and `2`.";
Options[GCTensorContract] = {PostETensorContract -> Identity};
GCTensorContract[t_, n1_List, n2_List, opt___] := Fold[GCTensorContract[#1, #2[[1]], #2[[2]], opt], t, Thread@{n1, n2}];
GCTensorContract[GCTensor[arr_, basis_], n1_Integer, n2_Integer, opt: OptionsPattern[]] := Module[
    {d, perm, trans},
    d = Length@basis;
    If[!UpDownChartsQ[#1, #2], Throw@Message[GCTensorContract::icpbs, #1, #2]] &[basis[[n1]], basis[[n2]]];
    perm = With[{range = Range@d}, InversePermutation@Join[Delete[range, {{n1}, {n2}}], range[[{n1, n2}]]]];
    DebugPrint["[GCTensorContrac] perm1 = ", perm];
    trans = GCTensorTranspose[GCTensor[arr, basis], perm];
    trans = EnsureGCTensorBasis[trans, Join[ConstantArray[None, Length@basis - 2], ChangeBasisIndex /@ basis[[{n1, n2}]]]];
    GCTensorContractLast2SameBasis[trans, OptionValue[PostETensorContract]]
];
SyntaxInformation[GCTensorContract] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

Options[GCTensorFixedContract] = Options@GCTensorContractTwo;
GCTensorFixedContract[GCTensor[arr_, basis_], metric_GCTensor, n_, opt: OptionsPattern[]] := GCTensorTranspose[
    GCTensorContractTwo[GCTensor[arr, basis], metric, {n}, {1}, Sequence@FilterRules[{opt}, Options@GCTensorContractTwo]],
    MoveTo[Length@basis, n]
];
SyntaxInformation[GCTensorFixedContract] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

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

ZeroGCTensorQ[GCTensor[expr_, {}]] := expr === 0;
ZeroGCTensorQ[GCTensor[expr_, {basis__}]] := With[{
    depth = Length@{basis},
    mapped = Map[ZeroETensorQ, expr, {Length@{basis}}]
}, And @@ If[depth > 1, Flatten[mapped, depth], mapped]];
SyntaxInformation[ZeroGCTensorQ] = {"ArgumentsPattern" -> {_}};

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
    ExcludeZeros -> True
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

DefGCTensorMapFunc[funcs__] := Function[func,
    GCTensor /: func[GCTensor[arr_, basis_], args___] := GCTensor[Map[func[#, args] &, arr, {Length@basis}], basis];
    GCTensor /: func[GCTensor[arr_, basis_][inds__], args___] := GCTensor[Map[func[#, args] &, arr, {Length@basis}], basis][inds];
] /@ {funcs};

DefGCTensorMapFunc[
    D,
    Dt,
    Simplify,
    ScreenDollarIndices,
    ToCanonical,
    NoScalar,
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

ValidateGCTensor[GCTensor[arr_, basis_]] := With[{
    subMs = SubManifoldsOfGChart@UpBasis@# & /@ basis,
    pms = If[MatchQ[#, _Symbol], 1, -1] & /@ basis
}];
SyntaxInformation[ValidateGCTensor] = {"ArgumentsPattern" -> {_}};

SignedVBundleOfIndex[id_Symbol] := VBundleOfIndex@id;
SignedVBundleOfIndex[-id_Symbol] := -VBundleOfIndex@id;

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
    }, If[Sort@ainds =!= Sort@ainds2, Throw@Message[CreateGCTensor::unmatched, ainds, ainds2]]];
    elem = If[Length@ainds > 0, ETensor[expr, ainds], expr];
    ReplacePart[arr, inds -> elem]
];
CreateGCTensor::unmatched = "Indices `1` and `2` don't match.";
CreateGCTensor[{components___}, charts_] := Module[
    {dim, arrDim, params, subVBundles, ret},
    dim = Length@charts;
    params = MapThread[If[MatchQ[#1, _Symbol], #2, -#2] &, {charts, CoordParamsOfGChart /@ charts}];
    subVBundles = SubVBundlesOfGChart /@ charts;
    ret = Fold[
        AddGCTensorElement[params, subVBundles],
        Array[InitGCTensorElement[Length /@ params, subVBundles], (Length /@ params) + (Length /@ subVBundles)],
        {components}
    ];
    GCTensor[ret, charts]
];
SyntaxInformation[CreateGCTensor] = {"ArgumentsPattern" -> {_, _}};

GCTensorPDGrad[chart_][expr_] := GCTensorPDGrad[expr, chart];
GCTensorPDGrad[expr_Plus, chart_] := GCTensorPDGrad[chart] /@ expr;
GCTensorPDGrad[GCTensor[arr_, basis_], chart_] := With[{
    ops = Join[
        ParamD[#[[3]]] & /@ CoordsOfGChart@chart,
        ETensorPDGrad@First@SlotsOfTensor@# & /@ SubManifoldsOfGChart@chart
    ],
    len = Length@basis
},
    GCTensor[Map[Function[elem, #@elem & /@ ops], arr, {len}], Append[basis, -chart]]
];
GCTensorPDGrad[scalar_, chart_] := With[{
    ops = Join[
        ParamD[#[[3]]] & /@ CoordsOfGChart@chart,
        ETensorPDGrad@First@SlotsOfTensor@# & /@ SubManifoldsOfGChart@chart
    ],
    len = Length@basis
},
    GCTensor[#@scalar & /@ ops, {-chart}]
];
SyntaxInformation[GCTensorPDGrad] = {"ArgumentsPattern" -> {_, _.}};

GCTensorPDDiv[GCTensor[arr_, basis_], n_Integer] := GCTensorPDDiv[GCTensor[arr, basis], n, basis[[n]]];
GCTensorPDDiv[GCTensor[arr_, basis_], n_Integer, basis2_Symbol] := With[{
    basis1 = basis[[n]],
    (* ckernels =  *)
}];

ContractTwoIndexedGCTensors[t1_GCTensor[inds1__], t2_GCTensor[inds2__], opt___] := Module[
    {pairs, pos1, pos2, res},
    pairs = xAct`xTensor`Private`TakePairs[{inds1}, {inds2}];
    pos1 = FirstPosition[{inds1}, #][[1]] & /@ pairs;
    pos2 = FirstPosition[{inds2}, ChangeIndex@#][[1]] & /@ pairs;
    res = GCTensorContractTwo[t1, t2, pos1, pos2, opt];
    res @@ Join[Delete[{inds1}, Transpose@{pos1}], Delete[{inds2}, Transpose@{pos2}]]
];

ContractOneIndexedGCTensors[t_GCTensor[inds__], opt___] := Module[
    {pairs, pos1, pos2, res},
    pairs = xAct`xTensor`Private`TakePairs@{inds};
    pos1 = FirstPosition[{inds}, #][[1]] & /@ pairs;
    pos2 = FirstPosition[{inds}, ChangeIndex@#][[1]] & /@ pairs;
    res = GCTensorContract[t, pos1, pos2, opt];
    res @@ Delete[{inds}, Join[Transpose@{pos1}, Transpose@{pos2}]]
];

ProductTwoIndexedGCTensors[t1_GCTensor[inds1__], t2_GCTensor[inds2]] := GCTensorProduct[t1, t2] @@ Join[{inds1}, {inds2}];

Null@{OneContraction, TwoContraction};
ShowContraction[OneContraction[_, idx_], tensors_] := "OneContraction"[List @@ tensors[[idx]]];
ShowContraction[TwoContraction[_, idx1_, idx2_], tensors_] := "TwoContraction"[List @@ tensors[[idx1]], List @@ tensors[[idx2]]];
ExecuteContraction[OneContraction[_, idx_], l_List, opt___] := {
    ContractOneIndexedGCTensors[l[[idx]], opt],
    Delete[l, idx]
};
ExecuteContraction[TwoContraction[_, idx1_, idx2_], l_List, opt___] := With[{
    t1 = l[[idx1]],
    t2 = l[[idx2]]
}, {
    ContractTwoIndexedGCTensors[t1, t2, opt],
    Delete[l, {{idx1}, {idx2}}]
}];
OptimizedGCTensorContraction[{}, opt___] = 1;
OptimizedGCTensorContraction[l_List, opt___] := Module[
    {tensors, contractions, contractionCounts, selected, res, l2},
    tensors = MapIndexed[If[MatchQ[#1, GCTensor[__][__]], {List @@ #1, #2[[1]]}, Nothing] &, l];
    If[Length@tensors > 0,
        contractions = Join[
            Table[If[Length@xAct`xTensor`Private`TakePairs[t[[1]], t[[1]]] > 0,
                OneContraction[
                    Length@xAct`xTensor`Private`DropPairs[t[[1]], t[[1]]],
                    t[[2]]
                ]
            ,
                Nothing
            ], {t, tensors}],
            Flatten[Table[TwoContraction[
                Length@xAct`xTensor`Private`DropPairs[tensors[[i, 1]], tensors[[j, 1]]],
                tensors[[i, 2]],
                tensors[[j, 2]]
            ], {i, 1, Length@tensors}, {j, i + 1, Length@tensors}], 2]
        ];
        If[Length@contractions > 0,
            contractionCounts = #[[1]] & /@ contractions;
            selected = contractions[[FirstPosition[contractionCounts, Min@contractionCounts][[1]]]];
            DebugPrint["selected: ", ShowContraction[selected, l]];
            {res, l2} = ExecuteContraction[selected, l, opt];
            OptimizedGCTensorContraction[Append[l2, res], opt]
        ,
            Times @@ l
        ]
    ,
        Times @@ l
    ]
];

Options[ContractGCTensors] = Union[Options[GCTensorContractTwo], Options[GCTensorContract]];
ContractGCTensors[covd_][expr_] := ContractGCTensors[expr, covd];
ContractGCTensors[expr_Plus, covd_, opt___] := ContractGCTensors[#, covd, opt] & /@ expr;
ContractGCTensors[expr_List, holder_, opt___] := ContractGCTensors[#, holder, opt] & /@ expr;
ContractNonTensors[e: GCTensor[__][__], __] := e;
ContractNonTensors[expr_, covd_, opt___] := ContractGCTensors[expr, covd, opt];
ContractGCTensors[expr_Times, covd_, opt___] := With[
    {l = ContractNonTensors[#, covd] & /@ List @@ expr},
    (Times @@ DeleteCases[l, GCTensor[__][__]]) * OptimizedGCTensorContraction[Cases[l, GCTensor[__][__]], opt]
];
ContractGCTensors[expr: GCTensor[__][__], _, opt___] := OptimizedGCTensorContraction[{expr}, opt];
ContractGCTensors[Power[a_, b_], holder_, opt___] := Power[ContractGCTensors[a, holder, opt], ContractGCTensors[b, holder, opt]];
ContractGCTensors[Scalar[expr_], holder_, opt___] := Scalar@ContractGCTensors[expr, holder, opt];
ContractGCTensors[expr_, __] := expr;
SyntaxInformation[ContractGCTensors] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

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

On[RuleDelayed::rhs];

End[];

Protect @@ Names[$Context <> "*"];

EndPackage[];