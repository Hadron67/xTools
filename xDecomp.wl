BeginPackage["xTools`xDecomp`", {"xAct`xCore`", "xAct`xTensor`", "xTools`xTension`"}];

(Unprotect[#]; ClearAll[#];) & /@ Names@{"`*", "`*"};

GChartQ::usage = "GChartQ[chart] gives True if chart is a generalized chart.";
DecompositionOfGChart::usage = "DecompositionOfGChart[chart] returns the form {n, {TangentM1, ...}}";
DecompositionOfSignedGChart::usage = "DecompositionOfSignedGChart[chart] is similar to DecompositionOfGChart[chart] but the subvbundles are signed accordingly.";
DimensionLabelsOfGChart::usage = "DimensionLabelsOfGChart[chart] returns the labels of each 1D submanifold of chart, used in e.g. CreateGCTensor.";
PDGChartInfo::usage = "PDGChartInfo[{p1, p2, ...}, {{CD1, TM1}, {CD2, TM2}, ...}, GCArray[...]] represents a PD object.";
PDInfoOfGChart::usage = "PDInfoOfGChart[chart] returns information for computing the partial derivative in the chart.";
BasisOfGChart::usage = "BasisOfGChart[chart] returns the basis vectors of the chart, in the form {{{et, edt}, ...}, {eX, ...}}.";
CoordBasisOfGChart::usage = "CoordBasisOfGChart[chart] gives single-variable coordinate list of the generalized chart, each element is in the form {e, ed, param}."
SubManifoldBasisOfGChart::usage = "SubManifoldBasisOfGChart[chart] gives a list of basis tensors of the generalized chart.";
MetricOfGChart::usage = "MetricOfGChart[chart] gives the metric of generalized chart, if defined."
GChartMetricDecompose::usage = "GChartMetricDecompose[chart, metric[-a, -b]] decomposes the metric using the defined metric decomposition rule.";
VBundleOfGChart::usage = "VBundleOfGChart[chart] gives the base VBundle of the generalized chart.";
SubVBundlesOfGChart::usage = "SubVBundlesOfGChart[chart] gives the subvbundles of sub manifolds of the generalized chart.";
IndexDimensionOfGChart::usage = "IndexDimensionOfGChart[chart] gives the index dimension, i.e., the dimension where sub manifolds are treated as dim 1, of the chart.";
GetBasisETensorsOfGChart::usage = "GetBasisETensorsOfGChart[chart] gives a list of basis of chart in ETensors form.";

DefGBasis::usage = "DefGBasis[chart, {{e, ed, param}...}, {eT...}] defines a generalized chart with coodinate basis and sub manifolds."
UndefGBasis::usage = "UndefGBasis[chart] undefines a generalized chart.";
DefGChart::usage = "DefGBasis[chart, vb, {coordCount, {Vvb1, vb2, ...}}] is a lightweight version of DefGBasis, where only coordinate and submanifold count are required.";
SetGChartMetric::usage = "SetGChartMetric[chart, metric, metricValue, metricInvValue] defines metric decomposition rule for metric.";

GetAllHeldTensors::usage = "GetAllCachedTensors[holder] gives a list of tensors that has a definition with the holder.";
GetAllHeldTensorRules::usage = "GetAllHeldTensorRules[holder] gives a rule list to replace all tensors defined in holder as the values.";
ClearGCTensorHolderCache::usage = "ClearGCTensorHolderCache[holder] clears the tensor cache.";

SetHeldMetric::usage = "SetHeldMetric[holder, metric, metric, invMetric] sets the metric of the holder.";
SetHeldPD::usage = "SetHeldPD[holder, chart] defines PD to holder.";
AddCurvatureTensorsToHolder::usage = "AddCurvatureTensorsToHolder[holder, chart, christoffelTensor] defines curvature tensors to the holder, using the given Christoffel tensor.";
SortGChartParamD::usage = "SortGChartParamD[expr, chart] or SortGChartParamD[chart][expr] makes all coordinate ParamDs of chart in expr insider PDs.";

HeldGCTensor::usage = "HeldGCTensor[holder, tensor] gives the tensor provider of the tensor.";
SlotsOfHeldGCTensor::usage = "SlotsOfHeldGCTensor[holder, tensor] gives the slots of the tensor returned by HeldGCTensor[holder, tensor], default to SlotsOfTensor[tensor].";
CachedGCTensor::usage = "CachedGCTensor[holder, tensor, inds] gives the cached tensor with specific risen/lowered indices.";
HeldMetricTensorQ::usage = "HeldMetricTensorQ[holder, metric] gives true if it holds the specified metric.";
HeldCovDOfGCTensorHolder::usage = "HeldCovDOfGCTensorHolder[holder, covd] returns the held CovD of covd.";

GCTensorHolderAction::usage = "GCTensorHolderAction[holder, tag, expr] defines some extra action applies in various computation stages. Typically used to simplify expressions.";
GCTensorHolderDefaultAction::usage = "GCTensorHolderDefaultAction[holder, tag, expr] is the default value for GCTensorHolderAction[holder, tag, expr].";
GCTensorHolderDAUseMetricVB::usage = "GCTensorHolderDAUseMetricVB[holder] is the UseMetricOnVBundle option passed to ToCanonical by GCTensorHolderDefaultAction[holder]. Default is All.";
GCTensorHolderDAContractMetric::usage = "GCTensorHolderDAContractMetric[holder] controls whether to run ContractMetric in GCTensorHolderDefaultAction[holder].";
PostCurvatureTensorCalculation::usage = "PostCurvatureTensorCalculation[tensor] is an action tag for GCTensorHolderAction that's executed after the calculation of every curvature related tensors.";
PostChangeIndex::usage = "PostChangeIndex[tensor, indices, changedInd] is an action tag for GCTensorHolderAction that's executed after changing index for a cached tensor.";

PDGChart::usage = "PDGChart[chart] represents the PD operator of the GChart chart.";
CovDGChart::usage = "CovDGChart[chart, chris] represents the CovD operator of the chart with Christoffel tensor chris.";
ExpandPDToBasis::usage = "ExpandPDToBasis[expr, chart] or ExpandPDToBasis[chart][expr] expands all partial derivatives in expr in terms of basis the coordinates and sub manifolds of chart.";
ExpandMetric::usage = "ExpandMetric[expr, chart] or ExpandMetric[chart][expr] replaces all metric in expr with the decomposed one defined in chart."
ExpandPDToGCTensor::usage = "ExpandPDToGCTensor[expr, chart] or ExpandPDToGCTensor[chart][expr] replaces all PDs to GCTensors in expr.";

GCTensor::usage = "GCTensor[values, charts] represents a generalized CTensor.";
GCArray::usage = "GCArray[arr, {{clen, subLen}, ...}] represents a generalized multi-dimension array.";
GCArrayBroadcastQ::usage = "GCArrayBroadcastQ[fn] gives true for functions that should be broadcasted into the GCTensor and GCArray.";
GCTensorChangeBasis::usage = "GCTensorChangeBasis[GCTensor[...], n, basis] changes the basis of n-th axis to the specified basis.";
EnsureGCTensorBasis::usage = "EnsureGCTensorBasis[GCTensor[...], basis] changes the basis of the GCTensor of any of the axes doesn't match the given one.";
GCTensorTranspose::usage = "GCTensorTranspose[GCTensor[...], perms] performs the tensor transpose of GCTensor.";
GCTensorContract::usage = "GCTensorContract[T, {{n1, n2}, ...}] contracts the n1 axis with n2 axis of T.";
GCTensorContractTwo::usage = "GCTensorContractTwo[T1, T2, {{n1, m1}, ...}] computes the tensor contraction between two tensors, in which n1-th axis of T1 and m1-th axis of T2 are contracted, and so on.";
GCTensorFixedContract::usage = "GCTensorFixedContract[T, T2, n] contracts n-th axis of T with the first axis of tensor T2, with the resulting axis stay at the original position. Typically used to change index.";
GCTensorToBasis::usage = "GCTensorToBasis[expr] converts GCTensor expressions to basis.";
GCTensorChangeIndices::usage = "GCTensorChangeIndices[GCTensor[...], inds, holder] changes the indices of GCTensor if unmatches with inds, using metric.";
ZeroGCTensorQ::usage = "ZeroGCTensorQ[T] gives True of GCTensor T is zero.";
PostETensorContract::usage = "PostETensorContract is an option of GCTensorContractTwo, GCTensorContract, ContractGCTensors that's called after each ETensor contraction.";
GCTensorToGCArray::usage = "GCTensorToGCArray[GCTensor[...]] converts the GCTensor into GCArray.";
GCArrayContract::usage = "GCArrayContract[T, {{n1, n2}, ...}] contracts the n1 axis with n2 axis of T";
GCArrayContractTwo::usage = "GCArrayContractTwo[T1, T2, {{n1, m1}, ...}] computes the tensor product of two arrays and then contracts n1 axes of T1 with m1 axes of T2.";
IdentityGCArray::usage = "IdentityGCArray[cLen, {TM1, ...}, {1, -1}] or IdentityGCArray[cLen, {TM1, ...}, {-1, 1}] returns an identity GCArray.";
GCArrayFromSparse::usage = "GCArrayFromSparse[{{...} -> ...}, {...}] initialize a GCArray from sparse notation.";
GCTensorFromSparse::usage = "GCTensorFromSparse[{{...} -> ...}, {chart1, ...}] initialize a GCTensor from sparse notation.";
RiemannOfPDGChartInfo::usage = "RiemannOfPDGChartInfo[PDGChartInfo[...]] gives the Riemann tensor of the connection defined by the PD.";
RiemannOfGChart::usage = "RiemannOfPDGChartInfo[chart] gives the GCTensor of RiemannOfPDGChartInfo[PDInfoOfGChart@chart].";

GCTensorToComponentRules::usage = "GCTensorToComponentRules[T] converts GCTensor to a list of component rules.";

ValidateGCTensor::usage = "ValidateGCTensor[GCTensor[...]] fixes some ";
ZeroGCTensor::usage = "ZeroGCTensor[charts] creats a GCTensor with zero components.";
CreateGCTensor::usage = "CreateGCTensor[{{r1, r2, ...} -> expr}, charts] is a convenient function for creating GCTensors.";
DeltaGCTensor::usage = "DeltaGCTensor[chart, {1, -1}] or DeltaGCTensor[chart, {-1, 1}] gives the delta tensor of the chart.";
GCTensorPDGrad::usage = "GCTensorPDGrad[expr, chart] calculates the partial derivative of the tensor expression in chart.";
GCTensorPDDiv::usage = "GCTensorPDDiv[expr, n, chart] acts the PD operator on expr and contracts it with n-th axis.";
GCTensorCovDGrad::usage = "GCTensorCovDGrad[expr, chart, chris] calculates the covariant derivative of the GCTensor using the given Christoffel tensor.";
GCTensorCovDDiv::usage = "GCTensorCovDDiv[expr, n, chart, chris] acts the CovD operator on expr and contracts it with n-th axis, using the given Christoffel tensor.";

ContractGCTensors::usage = "ContractGCTensors[expr, covd] performs contractions of GCTensor's, use metric if needed.";
ReplaceHeldGCTensors::usage = "ReplaceHeldGCTensors is an option of ContractGCTensors that controls whether to replace held tensors in the holder. Default is All.";
ReplaceHeldCovD::usage = "ReplaceHeldCovD is an option of ContractGCTensors that contains a list of CovDs to replace.";
OtherReplaces::usage = "OtherReplaces is a option of ContractGCTensors that defines additional tensors or covds to replace.";

Begin["`Private`"];

Off[RuleDelayed::rhs];

DecompositionOfGChart::undef = "GChart `1` is undefined.";
BasisOfGChart::undef = "GChart `1` has no defined basis.";
GChartQ[_] = False;
DecompositionOfGChart[a_] := (Message[DecompositionOfGChart::undef, a]; None);
DecompositionOfGChart[-a_?GChartQ] := DecompositionOfGChart[a];
DecompositionOfSignedGChart[a_?GChartQ] := DecompositionOfGChart[a];
DecompositionOfSignedGChart[-a_?GChartQ] := MapAt[-# &, DecompositionOfGChart[a], 2];
DimensionLabelsOfGChart[_] = {};
DimensionLabelsOfGChart[-a_?GChartQ] := DimensionLabelsOfGChart[a];
BasisOfGChart[a_] := (Message[BasisOfGChart::undef, a]; {});
BasisOfGChart[-a_?GChartQ] := BasisOfGChart[a];
PDInfoOfGChart[_] = None;
PDInfoOfGChart[-a_?GChartQ] := PDInfoOfGChart[a];
CoordsCountOfGChart[a_] := DecompositionOfGChart[a][[1]];
SubVBundlesOfGChart[chart_?GChartQ] := DecompositionOfGChart[chart][[2]];
SubVBundlesOfGChart[-chart_?GChartQ] := -SubVBundlesOfGChart[chart];
CoordBasisOfGChart[a_] := BasisOfGChart[a][[1]];
SubManifoldBasisOfGChart[a_] := BasisOfGChart[a][[2]];
VBundleOfGChart[-chart_?GChartQ] := VBundleOfGChart[chart];
SyntaxInformation[GChartQ] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[DecompositionOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[DecompositionOfSignedGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[CoordBasisOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SubManifoldBasisOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[VBundleOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SubVBundlesOfGChart] = {"ArgumentsPattern" -> {_}};

MetricOfGChart[_] = None;
GChartMetricDecompose[_, expr_] := expr;
SyntaxInformation[MetricOfGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[GChartMetricDecompose] = {"ArgumentsPattern" -> {_, _}};

IndexDimensionOfGChart[chart_?GChartQ] := DecompositionOfGChart[chart][[1]] + Length@DecompositionOfGChart[chart][[2]];
SyntaxInformation[IndexDimensionOfGChart] = {"ArgumentsPattern" -> {_}};

ChangeGBasisSign[a_?GChartQ] := -a;
ChangeGBasisSign[-a_?GChartQ] := a;
SignOfGBasis[a_?GChartQ] = 1;
SignOfGBasis[-a_?GChartQ] = -1;
DownGBasisQ[-_?GChartQ] = True;

DefGBasis[chart_Symbol, coordBasis_List, subManifolds_List] := Module[
    {params},
    GChartQ[chart] ^= True;
    VBundleOfGChart[chart] ^= VBundleOfIndex@subManifolds[[1, 2]];
    DimensionLabelsOfGChart[chart] ^= coordBasis[[All, 3]];
    DecompositionOfGChart[chart] ^= {Length@coordBasis, VBundleOfIndex@#[[1]] & /@ subManifolds};
    BasisOfGChart[chart] ^= {coordBasis[[All, {1, 2}, 0]], subManifolds[[All, 0]]};
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
    PDInfoOfGChart[chart] ^= PDGChartInfo[coordBasis[[All, 3]], Thread@{PD, SubVBundlesOfGChart@chart}];
];
SyntaxInformation[DefGBasis] = {"ArgumentsPattern" -> {_, _, _}};

UndefGBasis[l_List] := UndefGBasis /@ l;
UndefGBasis[chart_Symbol?GChartQ] := (
    Apply[{e, ed} |-> With[{
        Q = Symbol[ToString[First@SlotsOfTensor@e] <> "`Q"]
    },
        (* TODO: unset the following or remove it from definition *)
        (* param /: PD[-a_Symbol?Q][param] =.; *)
        UndefTensor /@ {e, ed};
    ], CoordBasisOfGChart@chart, {1}];
    UndefTensor /@ DecompositionOfGChart[chart][[2]];
    Remove@chart;
);
SyntaxInformation[UndefGBasis] = {"ArgumentsPattern" -> {___}};

DefGChart[chart_Symbol, vb_, {clen_, subs_}] := (
    GChartQ[chart] ^= True;
    VBundleOfGChart[chart] ^= vb;
    DecompositionOfGChart[chart] ^= {clen, subs};
    DimensionLabelsOfGChart[chart] ^= Range[clen + Length@subs];
);
DefGChart[chart_, vb_, li_, pd_PDGChartInfo] := (
    DefGChart[chart, vb, li];
    PDInfoOfGChart[chart] ^= pd;
);
SyntaxInformation[DefGChart] = {"ArgumentsPattern" -> {_, _, _, _.}};

PDGChartInfo[params_, subs_] := PDGChartInfo[params, subs, IdentityGCArray[Length@params, subs[[All, 2]], {-1, 1}]];
SyntaxInformation[PDGChartInfo] = {"ArgumentsPattern" -> {_, _.}};

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
    ETensor[#[[1]][ai], {ai}] & /@ CoordBasisOfGChart@chart,
    With[{bi = First@GetIndicesOfVBundle[First@SlotsOfTensor@#, 1]},
        ETensor[#[-bi, ai]]
    ] & /@ SubManifoldBasisOfGChart@chart
]];
GetBasisETensorsOfGChart[-chart_?GChartQ] := With[{
    ai = First@GetIndicesOfVBundle[VBundleOfGChart@chart, 1]
}, Join[
    ETensor[#[[2]][-ai], {-ai}] & /@ CoordBasisOfGChart@chart,
    With[{bi = First@GetIndicesOfVBundle[First@SlotsOfTensor@#, 1]},
        ETensor[#[bi, -ai]]
    ] & /@ SubManifoldBasisOfGChart@chart
]];
SyntaxInformation[GetBasisETensorsOfGChart] = {"ArgumentsPattern" -> {_}};

SignOfAIndex[xAct`xTensor`Labels] = 0;
SignOfAIndex[_LI] = 0;
SignOfAIndex[_Symbol] = 1;
SignOfAIndex[-_Symbol] = -1;
SignOfGChart[_?GChartQ] = 1;
SignOfGChart[-_?GChartQ] = -1;
UpGChart[a_?GChartQ] := a;
UpGChart[-a_?GChartQ] := a;

HeldGCTensor::undef = "tensor provider for `1` of holder `2` is not defined.";
HeldGCTensor[holder_, tensor_] := Null /; (Message[HeldGCTensor::undef, tensor, holder]; False);
SyntaxInformation[HeldGCTensor] = {"ArgumentsPattern" -> {_, _}};

SlotsOfHeldGCTensor[_, tensor_?xTensorQ] := SlotsOfTensor@tensor;
SyntaxInformation[SlotsOfHeldGCTensor] = {"ArgumentsPattern" -> {_, _}};

GCTensorHolderDAUseMetricVB[_] = All;
GCTensorHolderDAContractMetric[_] = True;
SyntaxInformation[GCTensorHolderDAUseMetricVB] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[GCTensorHolderDAContractMetric] = {"ArgumentsPattern" -> {_}};

runHolderContractMetric[holder_, expr_, args___] := If[GCTensorHolderDAContractMetric@holder, ContractMetric[expr, args], expr];
GCTensorHolderAction[holder_, tag_, expr_] := GCTensorHolderDefaultAction[holder, tag, expr];
GCTensorHolderDefaultAction[holder, _, e_] := e;
GCTensorHolderDefaultAction[holder_, PostChangeIndex[tensor_, inds_, i_], expr_] := Simplify@ToCanonical[runHolderContractMetric[holder, expr], UseMetricOnVBundle -> GCTensorHolderDAUseMetricVB@holder];
GCTensorHolderDefaultAction[holder_, PostChangeIndex[tensor_, inds_, i_], expr: GCTensor[arr_, basis_]] := Simplify@ToCanonical[
    runHolderContractMetric[holder, expr, First@MetricsOfVBundle@# & /@ SubVBundlesOfGChart@UpGChart@basis[[i]]],
    UseMetricOnVBundle -> GCTensorHolderDAUseMetricVB@holder
];
GCTensorHolderDefaultAction[holder_, PostCurvatureTensorCalculation[tensor_], expr_] := Simplify@ToCanonical[expr, UseMetricOnVBundle -> None];
SyntaxInformation[GCTensorHolderAction] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[GCTensorHolderDefaultAction] = {"ArgumentsPattern" -> {_, _, _}};

SyntaxInformation[PostChangeIndex] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[PostCurvatureTensorCalculation] = {"ArgumentsPattern" -> {_}};

CachedGCTensor[holder_, tensor_][inds___] := CachedGCTensor[holder, tensor, SignOfAIndex /@ {inds}][inds];
CachedGCTensor[holder_, tensor_, inds_] := (
    xToolsDebugPrint[CachedGCTensor, "calculating tensor ", tensor];
    holder /: CachedGCTensor[holder, tensor, inds] = HeldGCTensor[holder, tensor]
) /; inds === SignOfAIndex /@ SlotsOfHeldGCTensor[holder, tensor];
CachedGCTensor[holder_, tensor_, inds_] := With[{
    firstUnmatched = First@FirstPosition[MapThread[#1 =!= SignOfAIndex@#2 &, {inds, SlotsOfHeldGCTensor[holder, tensor]}], True],
    slots = SlotsOfHeldGCTensor[holder, tensor]
},
    holder /: CachedGCTensor[holder, tensor, inds] = With[{
        tensor2 = CachedGCTensor[holder, tensor, MapAt[-# &, inds, firstUnmatched]]
    },
        xToolsDebugPrint[CachedGCTensor, "calculating index change ", tensor, inds, " from ", MapAt[-# &, inds, firstUnmatched]];
        GCTensorFixedContract[
            tensor2,
            CachedGCTensor[holder, First@MetricsOfVBundle@slots[[firstUnmatched]], {#, #} &[inds[[firstUnmatched]]]],
            firstUnmatched
        ] // GCTensorHolderAction[holder, PostChangeIndex[tensor, inds, firstUnmatched], #] &
    ]
] /; inds =!= SignOfAIndex /@ SlotsOfHeldGCTensor[holder, tensor];
SyntaxInformation[CachedGCTensor] = {"ArgumentsPattern" -> {_, _, _.}};

MetricTensorOfGCTensorHolder[_] = None;
SyntaxInformation[MetricTensorOfGCTensorHolder] = {"ArgumentsPattern" -> {_}};

PersistentTensorCacheQ[_, _] = False;
HeldCovDOfGCTensorHolder[_, _] = None;
SetHeldMetric[holder_Symbol, metric_, metrici_GCTensor, invMetric_GCTensor] := (
    holder /: CachedGCTensor[holder, metric, {-1, -1}] = metrici;
    holder /: CachedGCTensor[holder, metric, {1, 1}] = invMetric;
    holder /: HeldMetricTensorQ[holder, metric] = True;
    holder /: PersistentTensorCacheQ[holder, metric] = True;
);
SetHeldPD[holder_Symbol, chart_?GChartQ] := (
    holder /: HeldCovDOfGCTensorHolder[holder, PD] = PDGChart[-chart];
);
SyntaxInformation[SetHeldMetric] = {"ArgumentsPattern" -> {_, _, _, _, _.}};
SyntaxInformation[SetHeldPD] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[HeldMetricTensorQ] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[HeldCovDOfGCTensorHolder] = {"ArgumentsPattern" -> {_, _}};

AddCurvatureTensorsToHolder[holder_Symbol, chart_?GChartQ, chris_] := Module[
    {vb, metric, cd, a0, b0, c0, d0},
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
            ContractGCTensors[ChristoffelToGradMetric[chris[a, -b, -c], metric] // ExpandPDToGCTensor[chart], holder]
        ,
            {a, -b, -c}
        ] // NoScalar
        // GCTensorHolderAction[holder, PostCurvatureTensorCalculation[chris], #] &;
        holder /: HeldCovDOfGCTensorHolder[holder, cd2] := CovDGChart[-chart, CachedGCTensor[holder, chris, {1, -1, -1}]];

        (* Riemann tensor *)
        holder /: HeldGCTensor[holder, riem] := RiemannOfGChart@chart + ETensor[
            ContractGCTensors[ChangeCurvature[riem[-a, -b, -c, d], cd2, PD] // ExpandPDToGCTensor[chart], holder]
        ,
            {-a, -b, -c, d}
        ] // NoScalar
        // SortCommParamDLeviCivitaCovD
        // ExpandParamDLeviCivitaChristoffel // ToCanonicalN
        // CovDCommuToRiemann // SeparateMetricRiemann // ToCanonicalN
        // GCTensorHolderAction[holder, PostCurvatureTensorCalculation[riem], #] &;
        holder /: SlotsOfHeldGCTensor[holder, riem] = MapAt[-# &, SlotsOfTensor@riem, 4];

        (* Ricci tensor *)
        holder /: HeldGCTensor[holder, ricci] := GCTensorContract[CachedGCTensor[holder, riem, {-1, -1, -1, 1}], {{2, 4}}] // ToCanonicalN
        // ContractTensorWithMetric[#, Riemann /@ subCDs, subMetrics] &
        // ToCanonicalN
        // GCTensorHolderAction[holder, PostCurvatureTensorCalculation[ricci], #] &;

        (* Ricci scalar *)
        holder /: HeldGCTensor[holder, ricciScalar] := GCTensor[
            ContractGCTensors[CachedGCTensor[holder, ricci][-a, -b] CachedGCTensor[holder, metric][a, b], holder]
            // ToCanonicalN
            // ContractTensorWithMetric[#, Ricci /@ subCDs, subMetrics] &
            // ToCanonicalN
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
    rules = (HoldPattern[ParamD[l___, #3, r___]@PD[-a_Symbol]@A_] :> ParamD[l, r]@PD[-a]@ParamD[#3]@A) & @@@ CoordBasisOfGChart@chart
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

ClearGCTensorHolderCache[holder_] := (Replace[{
    HoldPattern[Verbatim[HoldPattern][CachedGCTensor[holder, t_, inds_]] :> _] :> (holder /: CachedGCTensor[holder, t, inds] =.) /; !PersistentTensorCacheQ[holder, t]
}] /@ UpValues[holder]; Null);
SyntaxInformation[ClearGCTensorHolderCache] = {"ArgumentsPattern" -> {_}};

SyntaxInformation[PDGChart] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[CovDGChart] = {"ArgumentsPattern" -> {_, _}};

ExpandPDToBasis[chart_?GChartQ][expr_] := ExpandPDToBasis[expr, chart];
ExpandPDToBasis[expr_, chart_?GChartQ] := Module[
    {vb, a, A, replacement1, replacement2, rule},
    vb = VBundleOfGChart@chart;
    a = GetIndicesOfVBundle[vb, 1][[1]];
    replacement1[a_, A_] = Plus @@ (ParamD[#[[3]]][A] #[[2]][-a] & /@ CoordBasisOfGChart[chart]);
    replacement2[a_, A_] = With[{
        i1 = GetIndicesOfVBundle[First@SlotsOfTensor@#, 1][[1]],
        e1 = #
    }, {{i1}, e1[i1, -a] PD[-i1][A]}] & /@ SubManifoldBasisOfGChart[chart];
    rule = With[{
        Q = Symbol[ToString[vb] <> "`Q"]
    }, PD[-a_Symbol?Q][A_] :> replacement1[a, A] + Plus @@ Module @@@ replacement2[a, A]];
    expr /. rule
];
SyntaxInformation[ExpandPDToBasis] = {"ArgumentsPattern" -> {_, _.}};

ExpandPDToGCTensor[chart_][expr_] := ExpandPDToGCTensor[expr, chart];
ExpandPDToGCTensor[expr_, chart_?GChartQ] := With[{
    pmQ = Symbol[ToString[VBundleOfGChart@chart] <> "`pmQ"]
}, expr /. HoldPattern[PD[a_?pmQ]] :> PDGChart[-chart][a]];
SyntaxInformation[ExpandPDToGCTensor] = {"ArgumentsPattern" -> {_, _.}};

ExpandMetric[chart_?GChartQ][expr_] := ExpandMetric[expr, chart];
ExpandMetric[expr_, chart_?GChartQ] := With[{metric = MetricOfGChart@chart},
   expr /. HoldPattern[metric[inds__]] :> GChartMetricDecompose[chart, metric[inds]]
];
SyntaxInformation[ExpandMetric] = {"ArgumentsPattern" -> {_, _.}};

(* GCTensor *)
MoveTo[n_, n_] := {1};
MoveTo[from_Integer, to_Integer] := If[
    from < to,
    Join[Range[from - 1], {to}, Range[from, to - 1]],
    Join[Range[to - 1], Range[to + 1, from], {to}]
];
CompletePerm[perm_, len_] := With[{
    l = Length@perm
}, If[len > l, Join[perm, Range[l + 1, len]], perm]];

GCArrayTranspose[GCArray[array_, basis_], perms_] := With[{
    subMStartI = basis[[All, 1]],
    cperms = CompletePerm[perms, Length@basis]
},
    GCArray[
        If[cperms != Range@Length@basis,
            MapIndexed[Function[{elem, inds},
                ETensorTranspose[elem, InversePermutation@Ordering[
                    cperms[[
                        Position[inds - subMStartI, _?(# > 0 &), {1}][[All, 1]]
                    ]]
                ]]
            ], array, {Length@basis}] // Transpose[#, perms] &
        ,
            array
        ],
        Permute[basis, perms]
    ]
];

GCTensorTranspose[GCTensor[array_, basis_], perms_] := (
    GCArrayTranspose[GCTensorToGCArray@GCTensor[array, basis], perms][[1]] // GCTensor[#, Permute[basis, perms]] &
);
SyntaxInformation[GCTensorTranspose] = {"ArgumentsPattern" -> {_, _}};

GetBasisChangeMatrix[basisTo_, basisFrom_] := Outer[
    ETensorContractTwo[#1, #2, {{-1, -1}}] &,
    GetBasisETensorsOfGChart[basisTo],
    GetBasisETensorsOfGChart[ChangeGBasisSign@basisFrom]
];

Null@ContractPair;

GCTensorDotInner[] := GCTensorDotInner[Identity];
GCTensorDotInner[action_][lhs_, ContractPair[rhs_, len_]] := With[{r = Thread@{Range[-len, -1], Range[len]}},
    ETensorContractTwo[lhs, rhs, r] // If[len > 0, action@#, #] &
];

GCTensorChangeBasis::invalid = "Cannot change basis with different signs `1` to `2`.";
GCTensorChangeBasis[GCTensor[arr_, basis_], n_Integer, basis2_] := With[{
    delta1 = DeltaGCTensor[-basis[[n]], basis2]
},
    GCTensorFixedContract[GCTensor[arr, basis], DeltaGCTensor[-basis[[n]], basis2], n]
] /; If[SignOfGBasis[basis2] === SignOfGBasis@basis[[n]], True, Message[GCTensorChangeBasis::invalid, basis[[n]], basis2]; False];
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

UpDownGChartsQ[c_?GChartQ, -c_?GChartQ] = True;
UpDownGChartsQ[-c_?GChartQ, c_?GChartQ] = True;
UpDownGChartsQ[_, _] = False;

GCArrayDot[GCArray[arr1_, basis1_], GCArray[arr2_, basis2_], len_, action_] := Module[
    {coordLengths, a1, a2},
    coordLengths = basis1[[Length@basis1 - len + 1 ;;, 1]];
    a1 = If[len > 1, Map[Flatten, arr1, {Length[basis1] - len}], arr1];
    a2 = MapIndexed[Function[{elem, indices},
        Map[ContractPair[
            #,
            Length@Cases[Thread[indices > coordLengths], True]
        ] &, elem, {Length@basis2 - len}]
    ], arr2, {len}] // If[len > 1, Flatten[#, len - 1], #] &;
    GCArray[
        Inner[GCTensorDotInner[action], a1, a2, Plus],
        Join[
            basis1[[;; Length[basis1] - len]],
            basis2[[len + 1 ;;]]
        ]
    ]
];

GCArrayContractLast2[GCArray[arr_, basis_], len_, action_] := With[{
    ccount = basis[[-len ;; -1, 1]]
},
    Map[Function[mat,
        Total[MapIndexed[Function[{subMat, indices},
            With[{
                clen = Length@Cases[Thread[indices > ccount], True],
                elem = subMat[[Sequence @@ indices]]
            },
                ETensorContract[elem, Thread@{Range[-clen, -1], Range[-2 clen, -clen - 1]}] //
                    If[clen > 0, action@#, #] &
            ]
        ], mat, {len}], len]
    ], arr, {Length@basis - 2*len}] // GCArray[#, basis[[;; -2 len - 1]]] &
];

GChartListToGCArrayMeta[basis_] := With[{d = DecompositionOfGChart@#}, {d[[1]], Length@d[[2]]}] & /@ basis;
GCTensorToGCArray[GCTensor[arr_, basis_]] := GCArray[arr, GChartListToGCArrayMeta@basis];
SyntaxInformation[GCTensorToGCArray] = {"ArgumentsPattern" -> {_}};

Options[GCArrayContractTwo] = {PostETensorContract -> Identity};
GCArrayContractTwo[GCArray[arr1_, basis1_], GCArray[arr2_, basis2_], nn0_List, opt: OptionsPattern[]] := Module[
    {n1, n2, d1, d2, perm1, perm2, trans1, trans2},
    n1 = If[# < 0, Length@basis1 + # + 1, #] & /@ nn0[[All, 1]];
    n2 = If[# < 0, Length@basis2 + # + 1, #] & /@ nn0[[All, 2]];
    d1 = Length@basis1;
    d2 = Length@basis2;
    perm1 = InversePermutation@Join[Delete[Range@d1, Transpose@{n1}], n1];
    perm2 = InversePermutation@Join[n2, Delete[Range@d2, Transpose@{n2}]];
    xToolsDebugPrint[GCArrayContractTwo, "perm1 = ", perm1, ", perm2 = ", perm2];
    trans1 = GCArrayTranspose[GCArray[arr1, basis1], perm1];
    trans2 = GCArrayTranspose[GCArray[arr2, basis2], perm2];
    GCArrayDot[trans1, trans2, Length@n1, OptionValue[PostETensorContract]]
] /; Length@nn0 > 0;
GCArrayContractTwo[GCArray[arr1_, basis1_], GCArray[arr2_, basis2_], {}, opt: OptionsPattern[]] := GCArray[
    Map[Function[elem,
        Map[ETensorContractTwo[elem, #, {}] &, arr2, {Length@basis2}]
    ], arr1, {Length@basis1}],
    Join[basis1, basis2]
];
SyntaxInformation[GCArrayContractTwo] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

Options[GCArrayContract] = {PostETensorContract -> Identity};
GCArrayContract[GCArray[arr_, basis_], nn0_List, opt: OptionsPattern[]] := Module[
    {n1, n2, d, perm, trans},
    d = Length@basis;
    n1 = If[# < 0, Length@basis + # + 1, #] & /@ nn0[[All, 1]];
    n2 = If[# < 0, Length@basis + # + 1, #] & /@ nn0[[All, 2]];
    perm = With[{range = Range@d}, InversePermutation@Join[Delete[range, Transpose@{Join[n1, n2]}], range[[Join[n1, n2]]]]];
    xToolsDebugPrint[GCArrayContract, "perm1 = ", perm];
    trans = GCArrayTranspose[GCArray[arr, basis], perm];
    GCArrayContractLast2[GCArrayTranspose[GCArray[arr, basis], perm], Length@n1, OptionValue[PostETensorContract]]
] /; Length@nn0 > 0;
GCArrayContract[a_, {}, ___] := a;
SyntaxInformation[GCArrayContract] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

GCTensorContractTwo::icpbs = "cannot contract incompatible basis `1` and `2`.";
GCTensorContractTwo::argrx = "Expected 4 or more arguments, `1` found.";
Options[GCTensorContractTwo] = Join[Options[GCArrayContractTwo]];
GCTensorContractTwo[GCTensor[arr1_, basis1_], GCTensor[arr2_, basis2_], nn_List, opt: OptionsPattern[]] := With[{
    check = And @@ MapThread[If[!UpDownGChartsQ[#1, #2], Message[GCTensorContractTwo::icpbs, #1, #2]; False, True] &, {basis1[[ nn[[All, 1]] ]], basis2[[ nn[[All, 2]] ]]}]
},
    GCTensor[
        GCArrayContractTwo[
            GCTensorToGCArray@GCTensor[arr1, basis1],
            GCTensorToGCArray@EnsureGCTensorBasis[GCTensor[arr2, basis2], ReplacePart[
                ConstantArray[None, Length@basis2],
                #2 -> -basis1[[#1]] & @@@ nn
            ]],
            nn, opt
        ][[1]],
        Join[Delete[basis1, Transpose@{nn[[All, 1]]}], Delete[basis2, Transpose@{nn[[All, 2]]}]]
    ] /; check
];
GCTensorContractTwo[args___] := Null /; (
    xToolsDebugPrint[GCTensorContractTwo, "untransformed args: ", {args}];
    Message[GCTensorContractTwo::argrx, Length@{args}];
    False
);
SyntaxInformation[GCTensorContractTwo] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

GCTensorContract::icpbs = "cannot contract incompatible basis `1` and `2`.";
Options[GCTensorContract] = Options[GCArrayContract];
GCTensorContract[GCTensor[arr_, basis_], nn_List, opt: OptionsPattern[]] := With[{
    check = And @@ MapThread[If[!UpDownGChartsQ[#1, #2], Message[GCTensorContract::icpbs, #1, #2]; False, True] &, {basis[[ nn[[All, 1]] ]], basis[[ nn[[All, 2]] ]]}]
},
    GCTensor[
        GCArrayContract[
            GCTensorToGCArray@EnsureGCTensorBasis[
                GCTensor[arr, basis],
                ReplacePart[ConstantArray[None, Length@basis], #2 -> -basis[[#1]] & @@@ nn]
            ], nn, opt
        ][[1]],
        Delete[basis, Transpose@{Join @@ Transpose@nn}]
    ] /; check
];
SyntaxInformation[GCTensorContract] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

Options[GCTensorFixedContract] = Options@GCTensorContractTwo;
GCTensorFixedContract[GCTensor[arr_, basis_], metric_GCTensor, n_, opt: OptionsPattern[]] := GCTensorTranspose[
    GCTensorContractTwo[GCTensor[arr, basis], metric, {{n, 1}}, Sequence@FilterRules[{opt}, Options@GCTensorContractTwo]],
    MoveTo[Length@basis, n]
];
SyntaxInformation[GCTensorFixedContract] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

BasisETensorsFromSignedBasis[basis_Symbol] := Join[
    ToETensor@#[[1]] & /@ CoordBasisOfGChart@basis,
    ToETensor[#, {-1, 1}] & /@ SubManifoldBasisOfGChart@basis
];
BasisETensorsFromSignedBasis[-basis_Symbol] := Join[
    ToETensor@#[[2]] & /@ CoordBasisOfGChart@basis,
    ToETensor[#, {1, -1}] & /@ SubManifoldBasisOfGChart@basis
];
MoveIndicesOfVBundleToLeft[ETensor[expr_, inds_], vb_] := With[{
    pos = # & @@@ Position[vb === (VBundleOfIndex@#) & /@ inds, True, {1}]
}, ETensor[
    expr,
    Permute[inds, InversePermutation@Join[pos, Delete[Range@Length@inds, Transpose@{pos}]]]
]];
GCTensorToBasis[GCTensor[arr_, basis_]] := With[{
    bTensors = BasisETensorsFromSignedBasis /@ basis,
    clens = CoordsCountOfGChart@UpGChart@# & /@ basis,
    vb = VBundleOfGChart@UpGChart@First@basis
}, MapIndexed[Function[{elem, indices},
    ETensorContractTwo[
        MoveIndicesOfVBundleToLeft[ETensorProduct @@ MapThread[Part, {bTensors, indices}], vb],
        elem,
        Length@Position[Thread[indices > clens], True, {1}]
    ]
], arr, {Length@basis}] // Total[#, Length@basis] &];
SyntaxInformation[GCTensorToBasis] = {"ArgumentsPattern" -> {_}};

GCTensorChangeIndices[GCTensor[arr_, basis_], inds_, holder_] := Fold[
    GCTensorFixedContract[#1, #2[[3]], #2[[1]]] &,
    GCTensor[arr, basis],
    MapThread[
        With[{
            metric = First@MetricsOfVBundle@VBundleOfGChart@#3
        },
            If[HeldMetricTensorQ[holder, metric] && SignOfGChart[#3] != #2, {#1, #2, CachedGCTensor[holder, metric, {#2, #2}]}, Nothing]
        ] &
    , {Range@Length@basis, inds, basis}]
] /; Length@basis === Length@inds;
SyntaxInformation[GCTensorChangeIndices] = {"ArgumentsPattern" -> {_, _, _}};

ZeroGCTensorQ[GCTensor[expr_, {}]] := expr === 0;
ZeroGCTensorQ[GCTensor[expr_, {basis__}]] := With[{
    depth = Length@{basis},
    mapped = Map[ZeroETensorQ, expr, {Length@{basis}}]
}, And @@ If[depth > 1, Flatten[mapped, depth], mapped]];
SyntaxInformation[ZeroGCTensorQ] = {"ArgumentsPattern" -> {_}};

nonZeroComponentQ[_ -> ETensor[n_, __]] := n =!= 0;
nonZeroComponentQ[_ -> n_] := n =!= 0;
GCTensorToComponentRules[(GCTensor | GCArray)[arr_, basis_]] := Select[
    Flatten[MapIndexed[#2 -> CatchedScreenDollarIndices@#1 &, arr, {Length@basis}], Length@basis],
    nonZeroComponentQ
];
SyntaxInformation[GCTensorToComponentRules] = {"ArgumentsPattern" -> {_}};

GCArrayBroadcastQ@ScreenDollarIndices = True;
GCArrayBroadcastQ[fn_] := ETensorBroadcastQ@fn;
GCTensor /: fn_?GCArrayBroadcastQ[GCTensor[arr_, basis_], args___] := GCTensor[Map[fn[#, args] &, arr, {Length@basis}], basis];
GCTensor /: fn_?GCArrayBroadcastQ[GCTensor[arr_, basis_][inds__], args___] := GCTensor[Map[fn[#, args] &, arr, {Length@basis}], basis][inds];
GCArray /: fn_?GCArrayBroadcastQ[GCArray[arr_, basis_], args___] := GCArray[Map[fn[#, args] &, arr, {Length@basis}], basis];
GCArray /: fn_?GCArrayBroadcastQ[GCArray[arr_, basis_][inds__], args___] := GCArray[Map[fn[#, args] &, arr, {Length@basis}], basis][inds];

GCTensor[e_, {}][] := PutScalar@e;
GCTensor /: ETensor[t_GCTensor[inds__], {inds2__}] := With[{
    perm = PermutationProduct[InversePermutation@Ordering@{inds}, Ordering@{inds2}]
}, GCTensorTranspose[t, perm]] /; Sort@{inds} === Sort@{inds2};

GCTensor /: GCTensor[arr1_, basis_] + GCTensor[arr2_, basis_] := GCTensor[arr1 + arr2, basis];
GCTensor /: x_?xTools`xTension`Private`IndexedScalarQ * GCTensor[arr_, basis_] := GCTensor[arr * x, basis];
GCTensor /: GCTensor[arr1_, basis1_][inds1__] + GCTensor[arr2_, basis2_][inds2__] := With[{
    perm = PermutationProduct[InversePermutation@Ordering@{inds2}, Ordering@{inds1}]
},
    GCTensor[arr1 + GCArrayTranspose[GCTensorToGCArray@GCTensor[arr2, basis2], perm][[1]], basis1][inds1] /;
    Sort@{inds1} === Sort@{inds2} && basis1 === Permute[basis2, perm]
];

NonIndexedScalarQExt[ParamD[__]@expr_] := NonIndexedScalarQExt@expr;
NonIndexedScalarQExt[expr_] := xAct`xTensor`Private`NonIndexedScalarQ[expr];
GCTensor /: x_?NonIndexedScalarQExt * GCTensor[arr_, basis_][inds__] := GCTensor[arr * x, basis][inds];

SyntaxInformation[GCTensor] = {"ArgumentsPattern" -> {_, _}};

GCArray /: GCArray[arr1_, basis_] + GCArray[arr2_, basis_] := GCArray[arr1 + arr2, basis];
GCArray /: x_?xTools`xTension`Private`IndexedScalarQ * GCArray[arr_, basis_] := GCArray[arr * x, basis];

ValidateGCTensor[GCTensor[arr_, basis_]] := With[{
    subMs = SubManifoldBasisOfGChart@UpGChart@# & /@ basis,
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
        Length@p + FirstPosition[v, SignedVBundleOfIndex@s][[1]]
    ]]
, {spec, params, subVBundles}];
InitGCArrayElement[clens_, subVBundles_][indices__] := ZeroETensor@MapThread[PartOrNone, {subVBundles, {indices} - clens}];
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
    params = MapThread[If[MatchQ[#1, _Symbol], #2, -#2] &, {charts, DimensionLabelsOfGChart /@ charts}];
    subVBundles = SubVBundlesOfGChart /@ charts;
    ret = Fold[
        AddGCTensorElement[params, subVBundles],
        Array[InitGCArrayElement[Length /@ params, subVBundles], (Length /@ params) + (Length /@ subVBundles)],
        {components}
    ];
    GCTensor[ret, charts]
];
SyntaxInformation[CreateGCTensor] = {"ArgumentsPattern" -> {_, _}};

GCArrayFromSparse[elems_, basis_] := GCArray[
    ReplacePart[
        Array[InitGCArrayElement[basis[[All, 1]], basis[[All, 2]]], #1 + Length@#2 & @@@ basis],
        elems
    ], {#1, Length@#2} & @@@ basis
];
SyntaxInformation[GCArrayFromSparse] = {"ArgumentsPattern" -> {_, _}};

GCTensorFromSparse[elem_, basis_] := GCTensor[GCArrayFromSparse[elem, DecompositionOfSignedGChart /@ basis][[1]], basis];
SyntaxInformation[GCTensorFromSparse] = {"ArgumentsPattern" -> {_, _}};

pdChristOuter[{param_, n1_}, {{pd_, vb_}, n2_}] := With[{
    inds = MapAt[-# &, GetIndicesOfVBundle[vb, 3], {{2}, {3}}]
}, With[{
    elem = ETensor[-$RiemannSign * ParamD[param][Christoffel[pd] @@ inds], Permute[inds, {3, 1, 2}]]
}, {
    {n1, n2, n2, n2} -> elem,
    {n2, n1, n2, n2} -> -elem
}]];
RiemannOfPDGChartInfo[PDGChartInfo[coords_, subs_, mat_]] := Module[
    {coordsI, subsI, riems, pdChris},
    coordsI = MapIndexed[{#1, #2[[1]]} &, coords];
    subsI = MapIndexed[{#1, #2[[1]] + Length@coords} &, subs];
    riems = With[{
        idx = #2,
        pd = #1[[1]],
        inds = MapAt[-# &, -GetIndicesOfVBundle[#1[[2]], 4], 4]
    }, {idx, idx, idx, idx} -> ETensor[Riemann[pd] @@ inds, inds]] & @@@ subsI;
    pdChris = Flatten[Outer[pdChristOuter, coordsI, subsI, 1], 2];
    (* TODO: handle torsion *)
    GCArrayFromSparse[
        Join[riems, pdChris],
        With[{el = {Length@coords, subs[[All, 2]]}, el2 = {Length@coords, -subs[[All, 2]]}}, {el2, el2, el2, el}]
    ]
];
SyntaxInformation[RiemannOfPDGChartInfo] = {"ArgumentsPattern" -> {_}};

RiemannOfGChart[chart_] := GCTensor[RiemannOfPDGChartInfo[PDInfoOfGChart@chart][[1]], {-chart, -chart, -chart, chart}];
SyntaxInformation[RiemannOfGChart] = {"ArgumentsPattern" -> {_}};

IdentityGCArray[clen_, subvbs_, il: Alternatives[{1, -1}, {-1, 1}]] := With[{
    len = clen + Length@subvbs,
    l = {clen, Length@subvbs}
}, GCArray[
    ReplacePart[
        Array[InitGCArrayElement[{clen, clen}, il * {subvbs, subvbs}], {len, len}],
        Join[
            ({#, #} -> 1) & /@ Range@clen,
            MapIndexed[
                With[{
                    ai = il * GetIndicesOfVBundle[#, 2]
                },
                    {#2[[1]] + clen, #2[[1]] + clen} -> ETensor[delta @@ ai, ai]
                ] &
            ,
                subvbs
            ]
        ]
    ],
    {l, l}
]];
IdentityGCArray[clen_] := IdentityGCArray[clen, {}, {1, -1}];
SyntaxInformation[IdentityGCArray] = {"ArgumentsPattern" -> {_, _., _.}};

DeltaGCTensorSameBasis[chart_?GChartQ, il: Alternatives[{1, -1}, {-1, 1}]] := GCTensor[
    IdentityGCArray[CoordsCountOfGChart@chart, SubVBundlesOfGChart@chart, il][[1]],
    chart * il
];
MakeBasisChangeDeltaGCTensor[b1_, b2_] := GCTensor[Outer[
    ETensorContractTwo[#1, #2, {{-1, -1}}] &,
    GetBasisETensorsOfGChart[b1],
    GetBasisETensorsOfGChart[b2]
], {b1, b2}];
DeltaGCTensor[chart_?GChartQ][a_Symbol, -b_Symbol] := DeltaGCTensorSameBasis[chart, {1, -1}][a, -b];
DeltaGCTensor[chart_?GChartQ][-a_Symbol, b_Symbol] := DeltaGCTensorSameBasis[chart, {-1, 1}][-a, b];
DeltaGCTensor[chart_?GChartQ, -chart_?GChartQ] := DeltaGCTensorSameBasis[chart, {1, -1}];
DeltaGCTensor[-chart_?GChartQ, chart_?GChartQ] := DeltaGCTensorSameBasis[chart, {-1, 1}];
DeltaGCTensor[chart1_?GChartQ, -chart2_?GChartQ] := MakeBasisChangeDeltaGCTensor[chart1, -chart2];
DeltaGCTensor[-chart1_?GChartQ, chart2_?GChartQ] := MakeBasisChangeDeltaGCTensor[-chart1, chart2];
SyntaxInformation[DeltaGCTensor] = {"ArgumentsPattern" -> {_, _.}};

GCTensorPDGrad[chart_][expr_] := GCTensorPDGrad[expr, chart];
GCTensorPDGrad[expr_Plus, chart_] := GCTensorPDGrad[chart] /@ expr;
GCTensorPDGrad[expr_, -chart_?GChartQ] := GCTensorPDGrad[expr, -chart, PDInfoOfGChart@chart];
GCTensorPDGrad[GCTensor[arr_, basis_], -chart_?GChartQ, PDGChartInfo[coords_, subs_, mat_]] := With[{
    ops = Join[ParamD /@ coords, ETensorPDGrad @@@ subs]
}, GCTensor[
    GCArrayContractTwo[
        GCArray[
            Map[Function[elem, #@elem & /@ ops], arr, {Length@basis}],
            Append[GChartListToGCArrayMeta@basis, {Length@coords, Length@subs}]
        ],
        mat,
        {{-1, 2}}
    ][[1]],
    Append[basis, -chart]
]];
GCTensorPDGrad[scalar_, -chart_?GChartQ, PDGChartInfo[coords_, subs_, mat_]] := With[{
    ops = Join[ParamD /@ coords, ETensorPDGrad @@@ subs]
}, GCTensor[
    GCArrayContractTwo[
        mat,
        GCArray[#@scalar & /@ ops, {{Length@coords, Length@subs}}],
        {{2, 1}}
    ][[1]],
    {-chart}
]];
SyntaxInformation[GCTensorPDGrad] = {"ArgumentsPattern" -> {_, _.}};

GCTensorPDDiv[GCTensor[arr_, basis_], n_Integer] := GCTensorPDDiv[GCTensor[arr, basis], n, ChangeGBasisSign@basis[[n]]];
GCTensorPDDiv[GCTensor[arr_, basis_], n_Integer, -chart2_?GChartQ] := GCTensorContract[
    GCTensorPDGrad[GCTensor[arr, basis], -chart2],
    {{n, -1}}
];
SyntaxInformation[GCTensorPDDiv] = {"ArgumentsPattern" -> {_, _, _.}};

GCTensorCovDGrad[expr: GCTensor[arr_, basis_], -chart_?GChartQ, chris_] := With[{
    pd = GCTensorPDGrad[expr, -chart],
    len = Length@basis
},
    pd + Total@MapIndexed[With[{n = #2[[1]]},
        If[GChartQ[#],
            GCTensorTranspose[GCTensorContractTwo[expr, chris, {{n, 3}}], MoveTo[len, n]]
        ,
            GCTensorTranspose[GCTensorContractTwo[-expr, chris, {{n, 1}}], MoveTo[len + 1, n]]
        ]
    ] &, basis]
];
GCTensorCovDGrad[scalar_, -chart_?GChartQ, chris_] := GCTensorPDGrad[scalar, -chart]; (* TODO: torsion *)
SyntaxInformation[GCTensorCovDGrad] = {"ArgumentsPattern" -> {_, _, _}};

GCTensorCovDDiv[expr: GCTensor[arr_, basis_], n_, -chart_?GChartQ, chris_] := With[{
    pd = GCTensorPDDiv[expr, n, -chart],
    chris2 = GCTensorContract[chris, {{1, 2}}],
    len = Length@basis
},
    pd + Total@MapIndexed[With[{n2 = #2[[1]], cp = If[#2[[1]] > n, 1, 0]}, Which[
        n == n2,
        GCTensorContractTwo[expr, chris2, {{n2, 1}}],
        GChartQ[#1],
        GCTensorTranspose[
            GCTensorContractTwo[expr, chris, {{n2, 3}, {n, 2}}],
            MoveTo[len - 1, n2 - cp]
        ],
        True,
        GCTensorTranspose[
            GCTensorContractTwo[-expr, chris, {{n2, 1}, {n, 2}}],
            MoveTo[len - 1, n2 - cp]
        ]
    ]] &, basis]
];
SyntaxInformation[GCTensorCovDDiv] = {"ArgumentsPattern" -> {_, _, _, _}};

ContractTwoIndexedGCTensors[t1_GCTensor[inds1__], t2_GCTensor[inds2__], opt___] := Module[
    {pairs, pos1, pos2, res},
    pairs = xAct`xTensor`Private`TakePairs[{inds1}, {inds2}];
    pos1 = FirstPosition[{inds1}, #, None, {1}][[1]] & /@ pairs;
    pos2 = FirstPosition[{inds2}, ChangeIndex@#, None, {1}][[1]] & /@ pairs;
    res = GCTensorContractTwo[t1, t2, Thread@{pos1, pos2}, opt];
    res @@ Join[Delete[{inds1}, Transpose@{pos1}], Delete[{inds2}, Transpose@{pos2}]]
];

ContractOneIndexedGCTensors[t_GCTensor[inds__], opt___] := Module[
    {pairs, pos1, pos2, res},
    pairs = Union[UpIndex /@ xAct`xTensor`Private`TakePairs@{inds}];
    pos1 = FirstPosition[{inds}, #, None, {1}][[1]] & /@ pairs;
    pos2 = FirstPosition[{inds}, ChangeIndex@#, None, {1}][[1]] & /@ pairs;
    res = GCTensorContract[t, Thread@{pos1, pos2}, opt];
    res @@ Delete[{inds}, Join[Transpose@{pos1}, Transpose@{pos2}]]
];

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
            xToolsDebugPrint[OptimizedGCTensorContraction, "selected: ", ShowContraction[selected, l]];
            {res, l2} = ExecuteContraction[selected, l, opt];
            OptimizedGCTensorContraction[Append[l2, res], opt]
        ,
            Times @@ l
        ]
    ,
        Times @@ l
    ]
];

Options[GCTensorContractPrimitive] = {
    ReplaceHeldGCTensors -> All,
    ReplaceHeldCovD -> All,
    OtherReplaces -> {}
};
Options[ContractGCTensors] = Union[
    Options[GCTensorContractTwo],
    Options[GCTensorContract],
    Options[GCTensorContractPrimitive]
];
FilterSubContractOptions[opt___] := FilterRules[{opt}, Union[Options[GCTensorContractTwo], Options[GCTensorContract]]];
ContractGCTensors[covd_][expr_] := ContractGCTensors[expr, covd];
ContractGCTensors[expr_Plus, covd_, opt___] := ContractGCTensors[#, covd, opt] & /@ expr;
ContractGCTensors[expr_List, holder_, opt___] := ContractGCTensors[#, holder, opt] & /@ expr;
ContractGCTensors[expr_And, holder_, opt___] := ContractGCTensors[#, holder, opt] & /@ expr;
ContractGCTensors[expr_Or, holder_, opt___] := ContractGCTensors[#, holder, opt] & /@ expr;
ContractGCTensors[expr1_ == expr2_, holder_, opt___] := ContractGCTensors[expr1, holder, opt] == ContractGCTensors[expr1, holder, opt];
ContractGCTensors[Power[a_, b_], holder_, opt___] := Power[ContractGCTensors[a, holder, opt], ContractGCTensors[b, holder, opt]];
ContractGCTensors[Scalar[expr_], holder_, opt___] := Scalar@ContractGCTensors[expr, holder, opt];
ContractGCTensors[expr_, holder_, opt___] := OptimizedGCTensorContraction[{GCTensorContractPrimitive[expr, holder, Sequence@FilterRules[{opt}, Options[GCTensorContractPrimitive]]]}, Sequence@FilterSubContractOptions[opt]];

FilterReplaceGCTensor[tensor_, holder_, All] := MemberQ[GetAllHeldTensors[holder], tensor];
FilterReplaceGCTensor[tensor_, holder_, ts_List] := MemberQ[ts, tensor] && MemberQ[GetAllHeldTensors[holder], tensor];
FilterReplaceHeldCovD[covd_, holder_, All] := HeldCovDOfGCTensorHolder[holder, covd] =!= None;
FilterReplaceHeldCovD[covd_, holder_, td_List] := MemberQ[td, covd] && HeldCovDOfGCTensorHolder[holder, covd] =!= None;

GCTensorContractPrimitive[t_GCTensor[inds___], holder_, opt___] := GCTensorChangeIndices[t, SignOfAIndex /@ {inds}, holder][inds];
GCTensorContractPrimitive[PDGChart[a__][b__][expr__], holder_, opt___] := GCTensorEvalCovDGChart[PDGChart[a][b]@ContractGCTensors[expr, holder, opt], holder];
GCTensorContractPrimitive[CovDGChart[a1__][a2__][expr_], holder_, opt___] := GCTensorEvalCovDGChart[CovDGChart[a1][a2]@ContractGCTensors[expr, holder, opt], holder];
GCTensorContractPrimitive[t_?xTensorQ[inds___], holder_, opt: OptionsPattern[]] := GCTensorChangeIndices[t /. Select[OptionValue[OtherReplaces], #[[1]] === t &], SignOfAIndex /@ {inds}, holder][inds] /; FirstPosition[OptionValue[OtherReplaces], t -> _, None, {1}] =!= None;
GCTensorContractPrimitive[t_?xTensorQ[inds___], holder_, opt: OptionsPattern[]] := CachedGCTensor[holder, t][inds] /; FilterReplaceGCTensor[t, holder, OptionValue[ReplaceHeldGCTensors]];
GCTensorContractPrimitive[covd_?CovDQ[inds__][expr_], holder_, opt: OptionsPattern[]] := GCTensorContractPrimitive[
    HeldCovDOfGCTensorHolder[holder, covd][inds][expr],
    holder,
    opt
] /; FilterReplaceHeldCovD[covd, holder, OptionValue[ReplaceHeldCovD]];
GCTensorContractPrimitive[expr_, __] := expr;

GCTensorEvalCovDGChart[PDGChart[-chart_?GChartQ][-a_Symbol][t_GCTensor[l___, a_Symbol, r___]], holder_] := GCTensorPDDiv[t, Length@{l} + 1, -chart][l, r];
GCTensorEvalCovDGChart[PDGChart[-chart_?GChartQ][a_][t_GCTensor[inds__]], holder_] := GCTensorChangeIndices[GCTensorPDGrad[t, -chart], SignOfAIndex /@ {inds, a}, holder][inds, a] /; UpIndexQ[a] || !MemberQ[{inds}, ChangeIndex@a];
GCTensorEvalCovDGChart[CovDGChart[-chart_?GChartQ, chris_][-a_Symbol][t_GCTensor[l___, a_Symbol, r___]], holder_] := GCTensorCovDDiv[t, Length@{l} + 1, -chart, chris][l, r];
GCTensorEvalCovDGChart[CovDGChart[-chart_?GChartQ, chris_][a_][t_GCTensor[inds__]], holder_] := GCTensorChangeIndices[GCTensorCovDGrad[t, -chart, chris], SignOfAIndex /@ {inds, a}, holder][inds, a] /; UpIndexQ[a] || !MemberQ[{inds}, ChangeIndex@a];
GCTensorEvalCovDGChart[PDGChart[-chart_?GChartQ][a_][expr_?xAct`xTensor`Private`NonIndexedScalarQ], holder_] := GCTensorChangeIndices[GCTensorPDGrad[expr, -chart], {SignOfAIndex@a}, holder][a];
GCTensorEvalCovDGChart[CovDGChart[-chart_?GChartQ, chris_][a_][expr_?xAct`xTensor`Private`NonIndexedScalarQ], holder_] := GCTensorChangeIndices[GCTensorCovDGrad[expr, -chart, chris], {SignOfAIndex@a}, holder][a];

ContractGCTensors[expr_Times, covd_, opt___] := With[
    {l = ContractGCTensors[#, covd, opt] & /@ List @@ expr},
    (Times @@ DeleteCases[l, GCTensor[__][__]]) * OptimizedGCTensorContraction[Cases[l, GCTensor[__][__]], Sequence@FilterSubContractOptions[opt]]
];

SyntaxInformation[ContractGCTensors] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

On[RuleDelayed::rhs];

End[];

Protect @@ Names["`*"];

EndPackage[];