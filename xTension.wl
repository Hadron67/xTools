BeginPackage["xTools`xTension`", {"xAct`xCore`", "xAct`xTensor`", "xAct`xPerm`", "xAct`SymManipulator`"}];

(Unprotect[#]; ClearAll[#];) & /@ Names@{"`*", "`Private`*"};

ToCanonicalN::usage = "Calls ToCanonical with UseMetricOnVBundle -> None.";
SimplificationN::usage = "Calls Simplification[] using Implode.";

ErrorMarker::usage = "ErrorMarker[expr] is used to indicate an error input expr.";

CatchedScreenDollarIndices::usage = "CatchedScreenDollarIndices[expr] attempts to return ScreenDollarIndices[expr] and returns expr if failed.";
MakeDecomposedRules::usage = "MakeDecomposedRules[tensor[inds], values] generates rules by splitting all indices into subvbundles.";

SetDecomposedRules::usage = "SetDecomposedRules[tensor[inds], values] sets the generated decomposed rules to tensor.";

DefCoordinateParameter::usage = "DefCoordinateParameter[parentVB -> vb, param, e, ed] defines param as a coordinate parametre by defining e as the basis and ed as the dual basis, and define rules associated with them.";

BasisOfCoordinateParameter::usage = "BasisOfCoordinateParameter[param] returns the basis vector of coordinate parametre param.";
DualBasisOfCoordinateParameter::usage = "DualBasisOfCoordinateParameter[param] returns the basis dual vector of coordinate parametre param.";
CoordinateParameterOfVBundle::usage = "CoordinateParameterOfVBundle[vb] returns the coordinate parametre of vector bundle vb, where vb must be a subbundle of some bundle.";
VBundleOfCoordinateParameter::usage = "VBundleOfCoordinateParameter[param] returns the VBundle of coordinate parametre param";
CoordinateBasisQ::usage = "CoordinateBasisQ[e] returns True if e is defined as the coordinate basis or dual basis.";

DropCoordinateBasis::usage = "DropCoordinateBasis[expr] replaces all coordinate basis in expr with 1.";

SplitAllIndices::usage = "SplitAllIndices[expr, vb] or SplitAllIndices[vb][expr] splits all free indices in expr of VBundle vb.";

CalculateDecomposedChristoffel::usage = "CalculateDecomposedChristoffel[cd] calculates decomposed christoffel symbols of CovD cd in the form of rules.";

CalculateDecomposedRiemann::usage = "CalculateDecomposedRiemann[cd, chris] calculates decomposed riemann tensor using CovD cd and decomposed rules of Christoffel tensor chris.";

CalculateDecomposedRicci::usage = "CalculateDecomposedRicci[cd, riemannRules] calculates decomposed Ricci tensor using CovD cd and decomposed rules of Riemann tensor riemannRules.";

CalculateDecomposedRicciScalar::usage = "CalculateDecomposedRicciScalar[cd, ricciRules] calculates the decomposed Ricci scalar using CovD cd and decomposed Ricci tensor ricciRules.";

ChristoffelToRiemann::usage = "ChristoffelToRiemann[expr, cd] or ChristoffelToRiemann[cd][expr] tries to convert all Christoffel tensors in expr if CovD cd to Riemann tensor.";

DecomposedChristoffelRules::usage = "DecomposedChristoffelRules[cd] returns the decomposed Christoffel tensor rules of CovD cd, using cached result if already calculated.";

DecomposedRiemannRules::usage = "DecomposedRiemannRules[cd] returns the decomposed Riemann tensor rules of CovD cd, using cached result if already calculated.";

DecomposedRicciRules::usage = "DecomposedRicciRules[cd] returns the decomposed Ricci tensor rules of CovD cd, using cached result if already calculated.";

DecomposedRicciScalarRule::usage = "DecomposedRicciScalarRule[cd] returns the decomposed Ricci scalar rules of CovD cd, using cached result if already calculated.";

AllDecomposedRules::usage = "AllDecomposedRules[cd] returns all the decomposed rules of CovD cd, using cached result if already calculated.";

ClearDecomposedCache::usage = "ClearDecomposedCache[cd] clears all decomposed cache of Covd cd.";

ToDecomposed::usage = "ToDecomposed[expr, cd] or ToDecomposed[cd][expr] convert expr into decomposed forms of lists using rules AllDecomposedRules[cd].";

ReplaceIndicesRules::usage = "ReplaceIndicesRules[expr, fromVB, toVB, opt] returns a list of rules that replaces all indices of fromVB in expr to toVB.";
ReplaceCovDs::usage = "ReplaceCovDs[expr, cdFrom, cdTo, opt] replaces all indices and associated tensors (curvatures, etc) of cdFrom in expr with that of cdTo.";
DefRiemannVarD::usage = "DefRiemannVarD[cd] defines rules that makes VarD[Riemann[cd]][expr] work for expr containing Ricci tensors.";
UndefRiemannVarD::usage = "UndefRiemannVarD[cd] removes rules defined by DefRiemannVarD[cd]";
BHTemperature::usage = "BHTemperature[gtt, grr, gxx, r, r0] calculates the Hawking temperature of black hole with radial coordinate r and metric components gtt, grr, gxx, and horizon r0.";

SetCMetricRule::usage = "SetCMetricRule[metric, cmetric] defines component rules for metric using CTensor metric cmetric.";
WithCMetricRule::usage = "WithCMetricRule[expr_, metric_, cmetric_] evaluates expr by calling SetCMetricRule[metric, cmetric] first. It's also an option for ReplaceCovDs determining whether wrap the whole expression with WithCMetricRule or not.";
UnsetCMetricRule::usage = "UnsetCMetricRule[metric] removes all rules defined by SetCMetricRule[metric, ...].";
OtherRules::usage = "OtherRules is an option for ReplaceIndicesRules and ReplaceCovDs that defines other rules to be applied.";
MetricInv::usage = "MetricInv is an option for SetCMetricRule that specifies the inverse metric explicitly.";

DefMetricNsd::usage = "DefMetricNsd[metric[-a, -b], covd, ...] calls DefMetric[1, metric[-a, -b], covd, ...] and then defines SignDetOfMetric[metric] as poison.";
NoSignDet::usage = "NoSignDet can be used as the first argument of DefMetric to specify a metric with undefined signdet.";

ETensor::usage = "ETensor[expr, {freeIndices}] represents an expression with dummy and free indices.";
ETensorProduct::usage = "ETensorProduct[e1, e2, ...] computes the tensor outer product of ETensors e1, e2, ...";
ETensorTranspose::usage = "ETensorTranspose[ETensor[...], perms] performs the transpose of the ETensor.";
ETensorPDGrad::usage = "ETensorPDGrad[T, vb] or ETensorPDGrad[vb][T] adds a PD to the ETensor.";
ETensorPDDiv::usage = "ETensorPDDiv[T, axis] contracts the PD operator with the specfied axis.";

ZeroETensor::usage = "ZeroETensor[vbs] creates a zero ETensor.";
ETensorContract::usage = "ETensorContract[T, n1, n2] contracts n1 axis with n2 of T.";
ETensorContractTwo::usage = "ETensorContractTwo[T1, T2, {n1...}, {n2...}] contracts the n1... axes of T1 to n2.. axes of T2.";
ETensorRank::usage = "ETensorRank[ETensor[...]] gives the rank of the ETensor.";
ToETensor::usage = "ToETensor[T] converts tensor head T to ETensor.";
ZeroETensorQ::usage = "ZeroETensorQ[T] gives True of T is a zero tensor.";
UniqueIndex::usage = "UniqueIndex[a] gets a unique and temporary a-index of the form a$*.";
IndexRangeNS::usage = "IndexRangeNS[ns`a, ns`p] is similiar to IndexRange[a, p] except the returned indices are prefixed with namespace ns`.";
EulerDensityP::usage = "EulerDensityP[riem, D] gives the D dimension Euler density with Riemann tensor riem. Similar to xAct`xTras`EulerDensity except it does not multiply SigDet[metric].";
RiemannScalarList::usage = "RiemannScalarList[riem, n] gives a list of all possible scalars constructed from the Riemann tensor riem. RiemannScalarList[None, n] returns a list where Riemann tensors are represented by List.";

MakeCommParamDLeviCivitaCovD::usage = "MakeCommParamDLeviCivitaCovD[pd, cd[ind], expr] returns the commutator [pd, cd[ind]](expr) where pd is some derivative operator, assuming Levi-Civita connection.";
SortCommParamDLeviCivitaCovD::usage = "SortCommParamDLeviCivitaCovD[expr, filter] puts all ParamD inside CovD in expr, using commutator returned by MakeCommParamDLeviCivitaCovD.";
ExpandParamDLeviCivitaChristoffel::usage = "ExpandParamDLeviCivitaChristoffel[expr, filter] expands derivatives of the Christoffel tensor with Levi-Civita connection into derivatives of the metric.";
PdSymChristoffelToRiemann::usage = "PdSymChristoffelToRiemann[expr, filter] tries to recover Riemann tensors from torsion-free Christoffel tensors in expr.";
ChangeCovDNonChristoffel::usage = "ChangeCovDNonChristoffel[expr, cd1, cd2] applies ChangeCovD to tensors except Christoffel tensors.";
CovDCommuToRiemann::usage = "CovDCommuToRiemann[expr, filter] find covd commutators in expr and convert them into Riemann tensors.";
SeparateMetricRiemann::usage = "SeparateMetricRiemann[expr, filter] calls SeparateMetric on all Riemann tensors.";
ContractTensorWithMetric::usage = "ContractTensorWithMetric[expr, metric, filter] calls ContractMetric[t, metric] on specific tensors in expr.";

ListToCanonical::usage = "ListToCanonical[list, group] canonicalizes the list with the given symmetry group.";

$xToolsDebugFilter = {};
$xToolsDebugFilter::usage = "$xDecompDebugFilter is a global variable, containing all enabled debug messages.";

$xToolsDebugPrint = Print;
$xToolsDebugPrint::usage = "$xDecompDebugPrint is a global hook variable for printing debug messages.";
xToolsDebugPrint::usage = "xToolsDebugPrint[tag, msg] prints debug message with specified tag.";

Begin["`Private`"];

FilterExprList[All, __] = True;
FilterExprList[l_List, expr_, ___] := MemberQ[l, expr];
FilterExprList[func_, args__] := func[args];

xToolsDebugPrint[tag_, msg___] := If[FilterExprList[$xToolsDebugFilter, tag], $xToolsDebugPrint[tag, " ", msg]];
SetAttributes[DebugPrint, HoldRest];

CatchedScreenDollarIndices[expr_] := With[{v = Catch[ScreenDollarIndices@expr]}, If[v === Null, expr, v]];

(* metric decomposition *)

MakeDecomposedRules[tensor_[inds__], vals_] := Module[
    {spList, patList, ret, v, a, b, tensor2},
    spList = Tuples[Which[
        # === Labels, {#},
        VBundleQ[#] || VBundleQ[-#], List @@ First@SplittingsOfVBundle@UpIndex[#]
    ] & /@ SlotsOfTensor@tensor];
    patList = Function[it, Which[
        #2 === Labels, #1,
        VBundleQ[#2], #1 /. Verbatim[Pattern][a_, _] :> (PatternTest[b_, Evaluate@ToExpression[ToString@#2 <> "`Q"]] /. b->a)
    ] & @@@ Transpose@{{inds}, it}] /@ spList;
    ret = HoldPattern[#1] :> #2 & @@@ Transpose@{tensor2 @@@ patList, Flatten @ vals};
    ret /. tensor2 -> tensor
];
SyntaxInformation[MakeDecomposedRules] = {"ArgumentsPattern" -> {_, _}};

SetDecomposedRules[tensor_[inds__], vals_] := (#1 := #2) & @@@ MakeDecomposedRules[tensor[inds], vals];
SyntaxInformation[SetDecomposedRules] = {"ArgumentsPattern" -> {_, _}};

BasisOfCoordinateParameter[_] = Null;
DualBasisOfCoordinateParameter[_] = Null;
CoordinateParameterOfVBundle[_] = Null;
VBundleOfCoordinateParameter[_] = Null;
CoordinateBasisQ[_] = False;
SyntaxInformation[BasisOfCoordinateParameter] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[DualBasisOfCoordinateParameter] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[CoordinateParameterOfVBundle] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[VBundleOfCoordinateParameter] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[CoordinateBasisQ] = {"ArgumentsPattern" -> {_}};

DefCoordinateParameter::nsvb = "`1` is not a subvbundle of `2`";
DefCoordinateParameter[parentVB_?VBundleQ -> vb_?VBundleQ, param_?ParameterQ, e_, ed_] := Module[
    {a, b, c, otherVB},
    If[parentVB == vb || !SubvbundleQ[parentVB, vb], Throw[Message[DefCoordinateParameter::nsvb, vb, parentVB]]];

    {a} = GetIndicesOfVBundle[vb, 1];
    otherVB = DeleteCases[List @@ First@SplittingsOfVBundle@parentVB, vb];
    DefTensor[e[a], BaseOfVBundle@vb, PrintAs -> "(\!\(\*SubscriptBox[\(\[PartialD]\), \(" <> PrintAs[param] <> "\)]\))"];
    DefTensor[ed[-a], BaseOfVBundle@vb, PrintAs -> "(d" <> PrintAs[param] <> ")"];
    e /: e[b_Symbol] ed[-b_Symbol] = 1;
    CoordinateBasisQ[e] ^= True;
    CoordinateBasisQ[ed] ^= True;

    With[
        {pmQ = Symbol[ToString[#] <> "`pmQ"]},
        param /: PD[_?pmQ]@param = 0;
        e /: PD[_?pmQ]@e[_Symbol] = 0;
        ed /: PD[_?pmQ]@ed[-_Symbol] = 0;
    ] & /@ otherVB;
    With[
        {cd = CovDOfMetric@#, pmQ = Symbol[ToString[VBundleOfMetric@#] <> "`pmQ"]},
        param /: cd[_?pmQ]@param = 0;
        e /: cd[_?pmQ]@e[_Symbol] = 0;
        ed /: cd[_?pmQ]@ed[-_Symbol] = 0;
    ] & /@ Select[Union @@ MetricsOfVBundle /@ otherVB, CovDOfMetric[#] =!= PD &];

    CoordinateParameterOfVBundle[vb] ^= param;
    VBundleOfCoordinateParameter[param] ^= vb;
    BasisOfCoordinateParameter[param] ^= e;
    DualBasisOfCoordinateParameter[param] ^= ed;
];
SyntaxInformation[DefCoordinateParameter] = {"ArgumentsPattern" -> {_, _, _, _}};
DropCoordinateBasis[expr_] := expr /. _?CoordinateBasisQ[_] -> 1;
SyntaxInformation[DropCoordinateBasis] = {"ArgumentsPattern" -> {_}};

CoordinateParameterQ[p_] := BasisOfCoordinateParameter[p] =!= Null;
CoordinateParameterOfIndex[i_] := CoordinateParameterOfVBundle@VBundleOfIndex@UpIndex@i;
CoordinateIndexQ[i_] := AIndexQ[i] && CoordinateParameterOfIndex@i =!= Null;

Unprotect[PD, ParamD];
PD[-a_Symbol?CoordinateIndexQ][A_] := With[
    {param = CoordinateParameterOfIndex[a]},
    ParamD[param][A] DualBasisOfCoordinateParameter[param][-a]
];
ParamD[p_Symbol?CoordinateParameterQ]@PD[-a_Symbol]@A_ := PD[-a]@ParamD[p]@A;
Protect[PD, ParamD];

DecomposedIndicesList[ids_, vb_] := Module[
    {l, sp},
    l = Select[List @@ ids, ToExpression[ToString@vb <> "`pmQ"]];
    sp = List @@ First@SplittingsOfVBundle@vb;
    #1 -> IndexList @@ If[MatchQ[#1, _Symbol], #2, -#2] & @@@ Transpose@{l, Transpose[GetIndicesOfVBundle[#, Length@l] & /@ sp]}
];
ScalarTraceProductDummy[expr_] := expr /. Scalar[e_] :> Scalar@TraceProductDummy@e;
SplitAllIndices[l_List, vb_] := SplitAllIndices[#, vb] & /@ l;
SplitAllIndices[expr_, vb_] := expr // SplitIndex[#, DecomposedIndicesList[FindFreeIndices@#, vb]] & // ScalarTraceProductDummy // TraceProductDummy;
SplitAllIndices[vb_][expr_] := SplitAllIndices[expr, vb];
SyntaxInformation[SplitAllIndices] = {"ArgumentsPattern" -> {_, _.}};

ToCanonicalN[expr_, opt___] := ToCanonical[expr, UseMetricOnVBundle -> None, opt];
SimplificationN[expr_] := Implode[expr, ParamD] // Simplification // Explode;

CalculateDecomposed[expr_, vb_, calc_, calc2_] := Module[
    {freeInds, indsMap, vals, lhs},
    freeInds = FindFreeIndices@expr;
    indsMap = DecomposedIndicesList[freeInds, vb];
    vals = expr // calc // SplitIndex[#, indsMap] & // TraceProductDummy // Flatten // calc2 // SimplificationN;
    lhs = expr // SplitIndex[#, indsMap] & // Flatten;
    Flatten[MakeRule[#, MetricOn -> None] & /@ Transpose@{lhs, vals}]
];

SelectNonFlatCD[vbs_] := Select[CovDOfMetric /@ Select[Union @@ MetricsOfVBundle /@ vbs, !FlatMetricQ[#] &], v |-> v =!= Null];
CalculateDecomposedChristoffel[cd_] := Module[
    {metric, vb, sp, cds, a, b, c, vals},
    metric = MetricOfCovD@cd;
    vb = VBundleOfMetric@metric;
    sp = List @@ First@SplittingsOfVBundle@vb;
    cds = SelectNonFlatCD[sp];
    {a, b, c} = GetIndicesOfVBundle[vb, 3];
    CalculateDecomposed[
        Christoffel[cd][a, -b, -c],
        vb,
        ChristoffelToGradMetric,
        vals |-> Fold[ChangeCovD[#1, PD, #2] &, vals, cds]
    ]
];
SyntaxInformation[CalculateDecomposedChristoffel] = {"ArgumentsPattern" -> {_}};

ChristoffelToRiemann[cd_][expr_] := ChristoffelToRiemann[expr, cd];
ChristoffelToRiemann[expr_, cd_] := Module[
    {vb, a, b, c, d, chris, riemann, rule},
    {a, b, c, d} = GetIndicesOfVBundle[VBundleOfMetric@MetricOfCovD@cd, 4];
    chris = Christoffel@cd;
    riemann = Riemann@cd;
    rule = MakeRule[{
        Evaluate@PD[-b][chris[d, -a, -c]],
        Evaluate[-ChangeCurvature[riemann[-a, -b, -c, d], cd] +
        riemann[-a, -b, -c, d] + PD[-b][chris[d, -a, -c]]]
    }, MetricOn -> None, UseSymmetries -> False];
    ((expr /. rule) + expr) / 2
];
SyntaxInformation[ChristoffelToRiemann] = {"ArgumentsPattern" -> {_, _.}};

CalculateDecomposedRiemann[cd_, chrisRules_] := Module[
    {metric, vb, sp, cds, a, b, c, d},
    metric = MetricOfCovD@cd;
    vb = VBundleOfMetric@metric;
    sp = List @@ First@SplittingsOfVBundle@vb;
    {a, b, c, d} = GetIndicesOfVBundle[vb, 4];
    cds = SelectNonFlatCD@sp;
    CalculateDecomposed[
        Riemann[cd][-a, -b, -c, -d],
        vb,
        e |-> ChangeCurvature[e, cd],
        vals |-> Module[
            {v = vals /. chrisRules},
            v = Fold[ChristoffelToRiemann[#1, #2] &, v, cds];
            v = Fold[ChristoffelToGradMetric[#1, #2] &, v, MetricOfCovD /@ cds];
            Fold[ChangeCovD[#1, PD, #2] &, v, cds]
        ]
    ]
];
SyntaxInformation[CalculateDecomposedRiemann] = {"ArgumentsPattern" -> {_, _}};

CalculateDecomposedRicci[cd_, riemannRules_] := Module[
    {metric, riemann, ricci, a, b, c, d},
    metric = MetricOfCovD@cd;
    riemann = Riemann@cd;
    ricci = Ricci@cd;
    {a, b, c, d} = GetIndicesOfVBundle[VBundleOfMetric@metric, 4];
    CalculateDecomposed[ricci[-a, -b], VBundleOfMetric@metric, riemann[-a, -c, -b, -d] metric[c, d] &, e |-> e /. riemannRules]
];
SyntaxInformation[CalculateDecomposedRicci] = {"ArgumentsPattern" -> {_, _}};

CalculateDecomposedRicciScalar[cd_, ricciRules_] := Module[
    {metric, vb, a, b, ret},
    metric = MetricOfCovD@cd;
    vb = VBundleOfMetric@metric;
    {a, b} = GetIndicesOfVBundle[VBundleOfMetric@metric, 2];
    ret = With[
        {tmp = ATensor["tmp", -{vb, vb}]},
        tmp[-a, -b] metric[a, b] // TraceProductDummy // ReplaceAll[tmp -> Ricci[cd]] // ReplaceAll[ricciRules] // NoScalar // SimplificationN
    ];
    MakeRule[{Evaluate[RicciScalar[cd][]], Evaluate@Scalar@ret}, MetricOn -> None]
];
SyntaxInformation[CalculateDecomposedRicciScalar] = {"ArgumentsPattern" -> {_, _}};

DecomposedChristoffelRules[cd_] := DecomposedChristoffelRules[cd] ^= CalculateDecomposedChristoffel[cd];
DecomposedRiemannRules[cd_] := DecomposedRiemannRules[cd] ^= CalculateDecomposedRiemann[cd, DecomposedChristoffelRules[cd]];
DecomposedRicciRules[cd_] := DecomposedRicciRules[cd] ^= CalculateDecomposedRicci[cd, DecomposedRiemannRules[cd]];
DecomposedRicciScalarRule[cd_] := DecomposedRicciScalarRule[cd] ^= CalculateDecomposedRicciScalar[cd, DecomposedRicciRules[cd]];
SyntaxInformation[DecomposedChristoffelRules] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[DecomposedRiemannRules] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[DecomposedRicciRules] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[DecomposedRicciScalarRule] = {"ArgumentsPattern" -> {_}};
AllDecomposedRules[cd_] := Flatten@{
    DecomposedChristoffelRules[cd],
    DecomposedRiemannRules[cd],
    DecomposedRicciRules[cd],
    DecomposedRicciScalarRule[cd]
};
ClearDecomposedCache[cd_] := (
    cd /: DecomposedChristoffelRules[cd] =.;
    cd /: DecomposedRiemannRules[cd] =.;
    cd /: DecomposedRicciRules[cd] =.;
    cd /: DecomposedRicciScalarRule[cd] =.;
);
SyntaxInformation[AllDecomposedRules] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[ClearDecomposedCache] = {"ArgumentsPattern" -> {_}};

ToDecomposed[cd_][expr_] := ToDecomposed[expr, cd];
ToDecomposed[expr_List, cd_] := ToDecomposed[#, cd] & /@ expr;
ToDecomposed[expr_, cd_] := Module[
    {metric, sp, cds},
    metric = MetricOfCovD@cd;
    sp = List @@ First@SplittingsOfVBundle@VBundleOfMetric@metric;
    cds = Select[CovDOfMetric /@ Union @@ MetricsOfVBundle /@ sp, v |-> v =!= Null];
    expr // SeparateMetric[metric]
        // ChangeCovD[#, cd, PD] &
        // ReplaceDummies
        // SplitIndex[#, DecomposedIndicesList[FindFreeIndices@#, VBundleOfMetric@metric]] &
        // ScalarTraceProductDummy
        // TraceProductDummy
        // ReplaceAll@AllDecomposedRules[cd]
        // Fold[ChangeCovD[#1, PD, #2] &, #, cds] &
];
SyntaxInformation[ToDecomposed] = {"ArgumentsPattern" -> {_, _.}};

Options[ReplaceIndicesRules] = {OtherRules -> {}};
ReplaceIndicesRules::nsi = "index(ices) `1` not found in expression";
ReplaceIndicesRules[expr_, fromVB_, toVB_, opt: OptionsPattern[]] := Module[
    {inds, extraRules, notFoundInds},
    extraRules = OptionValue[OtherRules];
    inds = UpIndex /@ Select[List @@ FindIndices@expr, VBundleOfIndex[#] === fromVB &] // DeleteDuplicates;
    notFoundInds = DeleteCases[#1 & @@@ extraRules, Alternatives @@ inds];
    If[Length@notFoundInds > 0, Message[ReplaceIndicesRules::nsi, notFoundInds]];
    inds = DeleteCases[inds, Alternatives @@ (#1 & @@@ extraRules)];
    Join[extraRules, Thread[inds -> GetIndicesOfVBundle[toVB, Length@inds, #2 & @@@ extraRules]]]
];
SyntaxInformation[ReplaceIndicesRules] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

Options[ReplaceCovDs] = {OtherRules -> {}, WithCMetricRule -> False};
ReplaceCovDs[expr_, cdFrom_, cdTo_, opt: OptionsPattern[]] := With[{
    gFrom = MetricOfCovD@cdFrom,
    gTo = MetricOfCovD@cdTo,
    indRules = Select[OptionValue[OtherRules], !AIndexQ[#[[1]]] &],
    otherRules = Select[OptionValue[OtherRules], AIndexQ[#[[1]]] &],
    wrapper = If[OptionValue[WithCMetricRule],
        WithCMetricRule[#, First@MetricsOfVBundle@First@VBundlesOfCovD@cdTo, MetricOfCovD@cdTo] &,
        Identity
    ]
},
    wrapper@Unevaluated[expr /. Flatten@{
        indRules,
        ReplaceIndicesRules[
            expr,
            First@VBundlesOfCovD@cdFrom,
            First@VBundlesOfCovD@cdTo,
            OtherRules -> otherRules
        ],
        cdFrom -> cdTo,
        gFrom -> gTo,
        Riemann[cdFrom] -> Riemann[cdTo],
        Ricci[cdFrom] -> Ricci[cdTo],
        RicciScalar[cdFrom] -> RicciScalar[cdTo]
    }]
];
SyntaxInformation[ReplaceCovDs] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

DefRiemannVarD[cd_] := With[{
    metric = MetricOfCovD@cd,
    RiemannCD = Riemann@cd,
    RicciCD = Ricci@cd,
    RicciScalarCD = RicciScalar@cd,
    ii = GetIndicesOfVBundle[First@VBundlesOfCovD@cd, 2]
},
    RiemannCD /: ImplicitTensorDepQ[RicciCD, RiemannCD] = True;
    RiemannCD /: ImplicitTensorDepQ[RicciScalarCD, RiemannCD] = True;
    (RicciCD /: VarD[RiemannCD[inds__], cd2_][RicciCD[a_, b_], rest_] := Module[
        {#1, #2},
        VarD[RiemannCD[inds], cd2][RiemannCD[a, #1, b, #2], rest metric[-#1, -#2]]
    ]) & @@ ii;
    (RicciScalarCD /: VarD[RiemannCD[inds__], cd2_][RicciScalarCD[], rest_] := Module[
        {#1, #2},
        VarD[RiemannCD[inds], cd2][RicciCD[#1, #2], rest metric[-#1, -#2]]
    ]) & @@ ii;
];
SyntaxInformation[DefRiemannVarD] = {"ArgumentsPattern" -> {_}};
UndefRiemannVarD[cd_] := With[{
    RiemannCD = Riemann@cd,
    RicciCD = Ricci@cd,
    RicciScalarCD = RicciScalar@cd
},
    RiemannCD /: ImplicitTensorDepQ[RicciCD, RiemannCD] =.;
    RiemannCD /: ImplicitTensorDepQ[RicciScalarCD, RiemannCD] =.;
    RicciCD /: VarD[RiemannCD[inds__], cd2_][RicciCD[a_, b_], rest_] =.;
    RicciScalarCD /: VarD[RiemannCD[inds__], cd2_][RicciScalarCD[], rest_] =.;
];
SyntaxInformation[UndefRiemannVarD] = {"ArgumentsPattern" -> {_}};

Options[SetCMetricRule] = {MetricInv -> None};
SetCMetricRule[metric_, cmetric_, opt: OptionsPattern[]] := With[{
    inv = If[# =!= None, #, Inv@cmetric] &@OptionValue@MetricInv
},
    metric[a_?DownIndexQ, b_?DownIndexQ] := cmetric[a, b];
    metric[a_?UpIndexQ, b_?UpIndexQ] := inv[a, b];
];
SyntaxInformation[SetCMetricRule] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

UnsetCMetricRule[metric_] := (
    metric[a_?DownIndexQ, b_?DownIndexQ] =.;
    metric[a_?UpIndexQ, b_?UpIndexQ] =.;
);
SyntaxInformation[UnsetCMetricRule] = {"ArgumentsPattern" -> {_}};

Options[WithCMetricRule] = Options[SetCMetricRule];
WithCMetricRule[expr_, metric_, cmetric_, opt: OptionsPattern[]] := WithCleanup[
    SetCMetricRule[metric, cmetric, opt],
    expr,
    UnsetCMetricRule[metric]
];
SetAttributes[WithCMetricRule, HoldFirst];
SyntaxInformation[WithCMetricRule] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

BHTemperature[gtt_, grr_, gxx_, r_, r0_] := With[{fp = D[gtt/gxx, r], hp = D[1/grr, r]}, 1/(4 Pi) Sqrt[gxx fp hp] // Simplify // ReplaceAll[r -> r0]];
SyntaxInformation[BHTemperature] = {"ArgumentsPattern" -> {_, _, _, _, _}};

DefineTensorSyntaxInformation[T_[inds___], deps_, sym_, opt: OptionsPattern[]] := (
    SyntaxInformation[T] = {"ArgumentsPattern" -> Array[_ &, Length@{inds}]};
);
xTension["xTools`xTension`", DefTensor, "End"] := DefineTensorSyntaxInformation;

DefMetricNsd::nosigndet = "SignDet[`1`] is not defined.";
DefMetricNsd[metric_[inds___], covd_, args___] := (
    DefMetric[1, metric[inds], covd, args];
    SignDetOfMetric[metric] ^:= Throw@Message[DefMetricNsd::nosigndet, metric];
);
SyntaxInformation[DefMetricNsd] = {"ArgumentsPattern" -> {_, _, ___}};

NoSignDet::nsd = "SignDetOfMetric[`1`] is not defined.";
NoSignDet /: DefMetric[NoSignDet, metric_[inds___], args___] := (
    DefMetric[1, metric[inds], args];
    SignDetOfMetric[metric] ^:= Throw@Message[NoSignDet::nsd, metric];
);

SameIndexUDQ[a_Symbol, b_Symbol] = True;
SameIndexUDQ[-a_Symbol, -b_Symbol] = True;
SameIndexUDQ[_, _] = False;
CopyIndexSign[from_Symbol, to_Symbol] := to;
CopyIndexSign[-from_Symbol, to_Symbol] := -to;
CompatibleIndexListsQ[l1_, l2_] := Length@l1 === Length@l2 && Inner[SameIndexUDQ, l1, l2, And];

UniqueIndex[e_Symbol] := Module @@ ({{#}, #} &@xAct`xTensor`Private`NoDollar@e);
UniqueIndex[-e_Symbol] := -UniqueIndex[e];

(* ETensor *)
ETensor::invldmlt = "Attempting to multiply two ETensors.";
ETensor::icinds = "Incompatible indices `1` and `2`.";
ETensor[expr_] := ETensor[expr, List @@ FindFreeIndices@expr];
ETensor[expr_, {}] := Scalar[expr];
ETensor[expr_, inds_List][inds2__] /; Length@inds === Length@{inds2} := Module[
    {indPairs, matchedInds, unmatchedInds, unmatchedDummies, deltas},
    indPairs = Thread@{inds, {inds2}};
    matchedInds = Select[indPairs, SameIndexUDQ[#[[1]], #[[2]]] &];
    unmatchedInds = Append[#, UniqueIndex@#[[1]]] & /@ Select[indPairs, !SameIndexUDQ[#[[1]], #[[2]]] &];
    deltas = Times @@ (delta[ChangeIndex[#3], #2] & @@@ unmatchedInds);
    ReplaceIndex[ReplaceDummies@expr, Join[#1 -> #3 & @@@ unmatchedInds, #1 -> #2 & @@@ matchedInds]] * deltas
];
ETensor /:
    ETensor[expr1_, inds1_List] + ETensor[expr2_, inds2_List] := (
        ETensor[
            ReplaceDummies[expr1]
            + (ReplaceDummies[expr2] // ReplaceIndex[#, Thread[inds2 -> inds1]] &)
        , inds1]
    ) /; If[CompatibleIndexListsQ[inds1, inds2], True, Message[ETensor::icinds, inds1, inds2]; False];
ETensor /: Times[ETensor[expr_, inds_], factors__] := (
    If[Cases[{factors}, _ETensor], Message[ETensor::invldmlt]];
    ETensor[Times[expr, factors], inds]
);

DefETensorMapFunc[funcs__] := Function[func,
    ETensor /: func[ETensor[expr_, inds_], args___] := ETensor[func[expr, args], inds];
] /@ {funcs};
DefETensorMapFunc[
    D,
    Dt,
    ToCanonical,
    Simplification,
    Simplify,
    Factor,
    Together,
    ContractMetric,
    NoScalar,
    PdSymChristoffelToRiemann,
    SortCommParamDLeviCivitaCovD,
    ExpandParamDLeviCivitaChristoffel,
    ChangeCovDNonChristoffel,
    SeparateMetricRiemann,
    CovDCommuToRiemann,
    ContractTensorWithMetric
];
ETensor /: ParamD[params__][ETensor[expr_, args__]] := ETensor[ParamD[params][expr], args];
ETensor /: SeparateMetric[args___]@ETensor[expr_, args2__] := ETensor[SeparateMetric[args][expr], args2];
ETensor /: FindFreeIndices[ETensor[_, inds_]] := IndexList @@ inds;
ETensor /: ScreenDollarIndices[ETensor[expr_, inds_]] := Module[
    {dollars, rep},
    dollars = Select[Union[UpIndex /@ inds], MemberQ[Attributes@#, Temporary] &];
    rep = If[Length@dollars > 0, Thread[dollars -> GetIndicesOfVBundle[VBundleOfIndex@UpIndex@First@dollars, Length@dollars, Union[UpIndex /@ List @@ FindIndices@expr, UpIndex /@ inds]]], {}];
    ETensor[ScreenDollarIndices[expr /. rep], inds /. rep]
];
ETensor /: ChangeCovD[ETensor[expr_, inds_], args___] := ETensor[ChangeCovD[expr, args], inds];

SyntaxInformation[ETensor] = {"ArgumentsPattern" -> {_, _.}};

DedupeRules[inds1_, inds2_] := With[{
    dupes = Intersection[UpIndex /@ inds1, UpIndex /@ inds2]
}, With[{
    inds2d = Select[inds2, MemberQ[dupes, UpIndex@#] &]
}, Thread[inds2d -> (UniqueIndex /@ inds2d)]]];

IndexedScalarQ[_ETensor | _List] = False;
IndexedScalarQ[_] = True;

ETensorProduct[ETensor[expr1_, inds1_], ETensor[expr2_, inds2_], rest___] := With[{
    rep = DedupeRules[inds1, inds2]
},
    ETensorProduct[
        ETensor[ReplaceDummies[expr1] * ReplaceIndex[Evaluate@ReplaceDummies[expr2], rep], Join[inds1, inds2 /. rep]],
        rest
    ]
];
ETensorProduct[left___, Zero, right___] = Zero;
ETensorProduct[e_] := e;
ETensorProduct[expr_ /; !MatchQ[expr, _List], l_List, rest___] := ETensorProduct[expr, #, rest] & /@ l;
ETensorProduct[l_List, rest___] := ETensorProduct[#, rest] & /@ l;
ETensorProduct[ETensor[expr_, inds_], x_?IndexedScalarQ, rest___] := ETensorProduct[ETensor[expr * ReplaceDummies@x, inds], rest];
ETensorProduct[x_?IndexedScalarQ, ETensor[expr_, inds_], rest___] := ETensorProduct[ETensor[expr * ReplaceDummies@x, inds], rest];
ETensorProduct[x_?IndexedScalarQ, y_?IndexedScalarQ, rest___] := ReplaceDummies[x]*ReplaceDummies[y]*If[Length@{rest} === 0, 1, ETensorProduct[rest]];
SyntaxInformation[ETensorProduct] = {"ArgumentsPattern" -> {___}};

ETensorContractTwo0[expr1_, inds1_, expr2_, inds2_, n1_List, n2_List] := With[{
    a = inds1[[#]] & /@ n1,
    b = inds2[[#]] & /@ n2
},  With[{bi = UniqueIndex /@ b},
    ETensor[ReplaceDummies[expr1] * ReplaceIndex[ReplaceDummies[expr2], Thread[b -> bi]] * (Times @@ (delta[ChangeIndex@#1, ChangeIndex@#2] & @@@ Thread@{a, bi})), Join[Delete[inds1, Transpose@{n1}], Delete[inds2, Transpose@{n2}]]]
]];

ETensorContractTwo[ETensor[expr1_, inds1_], ETensor[expr2_, inds2_], n1_, n2_] := With[{
    rep = DedupeRules[inds1, inds2]
}, ETensorContractTwo0[expr1, inds1, ReplaceIndex[expr2, rep], inds2 /. rep, n1, n2]];
ETensorContractTwo[a_, b_, n_Integer] := ETensorContractTwo[a, b, Range[-n, -1], Range@n];
ETensorContractTwo[ETensor[expr_, inds_], x_?IndexedScalarQ, {}, {}] := ETensor[expr * ReplaceDummies@x, inds];
ETensorContractTwo[x_?IndexedScalarQ, ETensor[expr_, inds_], {}, {}] := ETensor[expr * ReplaceDummies@x, inds];
ETensorContractTwo[x_?IndexedScalarQ, y_?IndexedScalarQ, {}, {}] := ReplaceDummies[x] * ReplaceDummies[y];
SyntaxInformation[ETensorContractTwo] = {"ArgumentsPattern" -> {_, _, _, _}};

ETensorContract[ETensor[expr_, inds_], n1_List, n2_List] := With[{
    a1 = inds[[n1]],
    a2 = inds[[n2]]
}, With[{a2i = UniqueIndex /@ a2},
    ETensor[
        ReplaceIndex[expr, Thread[a2 -> a2i]] * (Times @@ MapThread[delta[ChangeIndex@#1, ChangeIndex@#2] &, {a1, a2i}]),
        Delete[inds, Transpose@{Join[n1, n2]}]
    ]
]];
ETensorContract[e_, {}, {}] := e;
SyntaxInformation[ETensorContract] = {"ArgumentsPattern" -> {_, _, _}};

ETensorRank[ETensor[_, inds_]] := Length@inds;
ETensorRank[_] = 0;
SyntaxInformation[ETensorRank] = {"ArgumentsPattern" -> {_}};

GetIndicesFromUpDownVBundle[vb_Symbol?VBundleQ, count_, exclude_] := GetIndicesOfVBundle[vb, count, exclude];
GetIndicesFromUpDownVBundle[-vb_Symbol?VBundleQ, count_, exclude_] := -GetIndicesOfVBundle[vb, count, exclude];
ToETensor[t_?xTensorQ] := With[{
    inds = Fold[Join[#1, GetIndicesFromUpDownVBundle[#2, 1, Union[UpIndex /@ #1]]] &, {}, SlotsOfTensor@t]
}, ETensor[t @@ inds, inds]];
ToETensor[t_?xTensorQ, indsSign_List] := With[{
    inds = Fold[Join[#1, GetIndicesOfVBundle[#2, 1, Union[UpIndex /@ #1]]] &, {}, SlotsOfTensor@t] * indsSign
}, ETensor[t @@ inds, inds]] /; Length@SlotsOfTensor@t === Length@indsSign;
SyntaxInformation[ToETensor] = {"ArgumentsPattern" -> {_, _.}};

ZeroETensorQ[ETensor[0, _]] = True;
ZeroETensorQ[0] = True;
ZeroETensorQ[_] = False;
SyntaxInformation[ZeroETensorQ] = {"ArgumentsPattern" -> {_}};

ETensorTranspose[ETensor[expr_, inds_], perms_] := ETensor[expr, Permute[inds, perms]];
ETensorTranspose[expr_, {}] := expr;
ETensorTranspose[expr_, {_Integer}] := expr;
ETensorTranspose[0, _] = 0;
ETensorTranspose[Zero, _] = Zero;
ETensorTranspose[expr_Plus, perms_] := ETensorTranspose[#, perms] & /@ expr;
SyntaxInformation[ETensorTranspose] = {"ArgumentsPattern" -> {_, _}};

ETensorPDGrad[vb_][expr_] := ETensorPDGrad[expr, vb];
ETensorPDGrad[x_?IndexedScalarQ, vb_] := With[{
    ai = First@GetIndicesOfVBundle[vb, 1]
}, ETensor[PD[-ai]@ReplaceDummies@x, {-ai}]];
ETensorPDGrad[ETensor[expr_, inds_], vb_] := With[{
    ai = UniqueIndex@First@GetIndicesOfVBundle[vb, 1]
}, ETensor[PD[-ai]@expr, Append[inds, -ai]]];
SyntaxInformation[ETensorPDGrad] = {"ArgumentsPattern" -> {_, _}};

ETensorPDDiv::iind = "Cannot contract index `1` with the derivative operator.";
ETensorPDDiv[n_][expr_] := ETensorPDDiv[expr, n];
ETensorPDDiv[ETensor[expr_, inds_], n_Integer] := With[{
    a1 = inds[[n]]
},
    If[!UpIndexQ@a1, Throw@Message[ETensorPDDiv::iind, ai]];
    ETensor[PD[-a1]@expr, Delete[inds, n]]
];
SyntaxInformation[ETensorPDDiv] = {"ArgumentsPattern" -> {_, _}};

ZeroETensor[{}] = 0;
ZeroETensor[vbs_] := ETensor[0, Fold[With[{
    ind = GetIndicesOfVBundle[#2, 1, UpIndex /@ #1]
}, Join[#1, If[MatchQ[#2, _Symbol], ind, -ind]]] &, {}, vbs]];
SyntaxInformation[ZeroETensor] = {"ArgumentsPattern" -> {_}};

IndexRangeNS::nsne = "Symbols have different namespaces `1` and `2`.";
IndexRangeNS[a_Symbol, p_Symbol] := With[{
    ns = Context@a,
    a0 = SymbolName@a,
    p0 = SymbolName@p
}, Symbol[ns <> #] & /@ CharacterRange[a0, p0]] /; Context@a === Context@p;
IndexRangeNS[a_Symbol, p_Symbol] := Throw@Message[IndexRangeNS::nsne, Context@a, Context@p] /; Context@a =!= Context@p;
SyntaxInformation[IndexRangeNS] = {"ArgumentsPattern" -> {_, _}};

EulerDensityP[riem_, dim_?EvenQ] := With[{
    inds = GetIndicesOfVBundle[-First@SlotsOfTensor@riem, 2dim],
    rangDim = Range@dim
}, With[{
    riemProd = Product[riem[Slot[2i - 1], Slot[2i], -inds[[2i - 1]], -inds[[2i]]], {i, 1, dim / 2}]
},
    SameDummies@Total[Signature@# * ToCanonical[riemProd & @@#] &@PermuteList[inds, InversePerm@#] & /@ TransversalInSymmetricGroup[
        StrongGenSet[rangDim, GenSet @@ (Cycles /@ Partition[rangDim, 2])],
        Symmetric@rangDim
    ]]
]];
SyntaxInformation[EulerDensityP] = {"ArgumentsPattern" -> {_, _}};

(* RiemannScalarList *)

RiemTermCanonicalQ[{None, None, None, None}] = True;
RiemTermCanonicalQ[{_Integer, None, None, None}] = True;
RiemTermCanonicalQ[{_Integer, None, _Integer, None}] = True;
RiemTermCanonicalQ[{a_Integer, b_Integer, None, None}] := a != -b;
RiemTermCanonicalQ[{a_Integer, b_Integer, c_Integer, None}] := a != -b && (a < 0 || b != -c);
RiemTermCanonicalQ[{a_Integer, b_Integer, c_Integer, d_Integer}] := a != -b && c != -d;
RiemTermCanonicalQ[_] = False;
AtSlot12Or34Q[{n_Integer, 1}, {n_Integer, 2}] = True;
AtSlot12Or34Q[{n_Integer, 2}, {n_Integer, 1}] = True;
AtSlot12Or34Q[{n_Integer, 3}, {n_Integer, 4}] = True;
AtSlot12Or34Q[{n_Integer, 4}, {n_Integer, 3}] = True;
AtSlot12Or34Q[__] = False;
ContainsMixedContractTermQ[{a_, b_, c_, d_}, list_] := (
    (a != None && c != None && AtSlot12Or34Q[
        FirstPosition[list, -a, None, {2}][[2]],
        FirstPosition[list, -c, None, {2}][[2]]
    ]) && (b != None && d != None && AtSlot12Or34Q[
        FirstPosition[list, -b, None, {2}][[2]],
        FirstPosition[list, -d, None, {2}][[2]]
    ])
);
ContainsMixedContractTermQ[{a_, b_, c_, d_}, list_] := With[{
    apos = FirstPosition[list, -a, None, {2}],
    bpos = FirstPosition[list, -b, None, {2}],
    cpos = FirstPosition[list, -c, None, {2}],
    dpos = FirstPosition[list, -d, None, {2}]
},
    (a =!= None && c =!= None && AtSlot12Or34Q[apos, cpos])
    || (b =!= None && d =!= None && AtSlot12Or34Q[bpos, dpos])
    || (a =!= None && d =!= None && AtSlot12Or34Q[apos, dpos])
    || (b =!= None && c =!= None && AtSlot12Or34Q[bpos, cpos])
];
MixedContractTermQ[list_] := Or @@ (ContainsMixedContractTermQ[#, list] & /@ list);
RiemTermOuterLineCount[l_List] := With[{
    l2 = DeleteCases[l, None]
},
    Length@Complement[l2, -l2]
];
RiemCanonicalQ[l_List] := And[
    And @@ (RiemTermCanonicalQ /@ l),
    OrderedQ[RiemTermOuterLineCount /@ l, GreaterEqual],
    !MixedContractTermQ[l]
];
NextContractNum[{}] = 1;
NextContractNum[l_List] := Max@l + 1;
RiemAddContraction[list_] := With[{
    n = NextContractNum@DeleteCases[Flatten@list, None],
    np = Position[list, None, {2}]
},
    Select[
        ReplacePart[list, {np[[1]] -> n, # -> -n}] & /@ np[[2 ;;]],
        RiemCanonicalQ
    ]
];
ActiveRiem[l_] := l /; FreeQ[l, None];
RiemEnumerateContractions[n_] := {ActiveRiem[ConstantArray[None, {n, 4}]]} //. HoldPattern[ActiveRiem[l_]] :> (Sequence @@ (ActiveRiem /@ RiemAddContraction[l]));

RiemannScalarList[None, n_] := RiemEnumerateContractions[n];
SyntaxInformation[RiemannScalarList] = {"ArgumentsPattern" -> {_, _}};

(* Christoffel derivatives *)

PdLeviCivitaGamma[pd_, cd_, metric_, a_?UpIndexQ, -b_?UpIndexQ, -c_?UpIndexQ] := With[
    {d = UniqueIndex@a},
    1/2 metric[a, d] (cd[-b]@pd@metric[-d, -c] + cd[-c]@pd@metric[-b, -d] - cd[-d]@pd@metric[-b, -c])
];

MakeCommParamDLeviCivitaCovD[pd_, cd_[-cdInd_?UpIndexQ], expr_] := With[{
    metric = MetricOfCovD@cd,
    i = UniqueIndex@cdInd,
    cdvb = VBundleOfIndex@cdInd
},
    Total[If[UpIndexQ@#,
       PdLeviCivitaGamma[pd, cd, metric, #, -cdInd, -i] ReplaceIndex[expr, {# -> i}]
    ,
       -PdLeviCivitaGamma[pd, cd, metric, i, -cdInd, #] ReplaceIndex[expr, {# -> -i}]
    ] & /@ Select[List @@ FindFreeIndices@expr, VBundleOfIndex@# === cdvb &]]
];
MakeCommParamDLeviCivitaCovD[pd_, cd_[cdInd_?UpIndexQ], expr_] := With[{
    metric = MetricOfCovD@cd,
    i = UniqueIndex@cdInd,
    i2 = UniqueIndex@cdInd
},
    metric[cdInd, i] MakeCommParamDLeviCivitaCovD[pd, cd[-i], expr] - metric[cdInd, i] pd@metric[-i, -i2] cd[i2]@expr
];
SyntaxInformation[MakeCommParamDLeviCivitaCovD] = {"ArgumentsPattern" -> {_, _, _}};

SortCommParamDLeviCivitaCovD[expr_, filter_] := expr //. {
    ParamD[pl___, p_]@cd_?CovDQ[cdInd_]@e2_ :> ParamD[pl][
        cd[cdInd]@ParamD[p]@e2 + MakeCommParamDLeviCivitaCovD[ParamD[p], cd[cdInd], e2]
    ] /; FilterExprList[filter, cd]
};
SortCommParamDLeviCivitaCovD[expr_] := SortCommParamDLeviCivitaCovD[expr, All];
SyntaxInformation[SortCommParamDLeviCivitaCovD] = {"ArgumentsPattern" -> {_, _.}};

ChristoffelPDQ[t_?xTensorQ] := ContainsAll[TensorID@t, {Christoffel, PD}];
ChristoffelPDQ[_] = False;
SymChristoffelQ[t_?xTensorQ] := MemberQ[TensorID@t, Christoffel] && SymmetryGroupOfTensor@t === StrongGenSet[{2, 3}, GenSet@xAct`xPerm`Cycles@{2, 3}];
SymChristoffelQ[_] = False;
NonChristoffelQ[t_?xTensorQ] := !MemberQ[TensorID@t, Christoffel];
NonChristoffelQ[_] = True;
CovDOfChristoffelPD[t_] := First@DeleteCases[Cases[TensorID@t, _?CovDQ], PD];
RiemannTensorQ[t_?xTensorQ] := MemberQ[TensorID@t, Riemann];
RiemannTensorQ[_] = False;
RicciTensorQ[t_?xTensorQ] := MemberQ[TensorID@t, Ricci];
RicciTensorQ[_] = False;

ExpandParamDLeviCivitaChristoffel[expr_, filter_] := expr /. HoldPattern[
    ParamD[pl___, p_]@chris_?ChristoffelPDQ[a_?UpIndexQ, -b_?UpIndexQ, -c_?UpIndexQ]
] :> With[{
        cd = CovDOfChristoffelPD@chris
    }, ParamD[pl]@PdLeviCivitaGamma[ParamD[p], cd, MetricOfCovD@cd, a, -b, -c]
] /; FilterExprList[filter, chris];
ExpandParamDLeviCivitaChristoffel[expr_] := ExpandParamDLeviCivitaChristoffel[expr, All];
SyntaxInformation[ExpandParamDLeviCivitaChristoffel] = {"ArgumentsPattern" -> {_, _.}};

UpDownRules[l_List] := Replace[a_ -> b_] /@ l;

TestPdSymChrisPair[exprL_List][{n1_, chris_, {a1_, b1_, c1_, d1_}, _}, {n2_, chris_, _, _}] := With[{
    expr1 = exprL[[n1, 2]],
    expr2 = exprL[[n2, 2]]
}, If[n1 == n2, None, With[{
    sym1 = ReplaceIndex[expr1, {a1 -> b1, -a1 -> -b1, b1 -> a1, -b1 -> -a1}],
    sym2 = ReplaceIndex[expr1, {a1 -> c1, -a1 -> -c1, c1 -> a1, -c1 -> -a1}]
}, Which[
    ToCanonical[expr2 + sym1, UseMetricOnVBundle -> None] === 0, 1,
    ToCanonical[expr2 + sym2, UseMetricOnVBundle -> None] === 0, 2,
    True, None
]]]];
PdSymChrisPairToSymChris2[{_, chris_, {a1_, b1_, c1_, d1_}, factor_}, w_] := With[{
    e1 = UniqueIndex@a1
}, With[{
    remnant = chris[d1, -b1, -e1] chris[e1, -a1, -c1] - chris[d1, -a1, -e1] chris[e1, -b1, -c1] - Riemann[CovDOfChristoffelPD@chris][-a1, -b1, -c1, d1]
}, If[w == 1, remnant, remnant /. {b1 -> c1, c1 -> b1}] * factor]];
FoldBinaryTest[operator_, list_] := Fold[Function[{ret, pair},
    With[{
        elem1 = list[[pair[[1]]]],
        elem2 = list[[pair[[2]]]],
        selectedSet = Union @@ ret[[All, 1 ;; 2]]
    }, If[MemberQ[selectedSet, elem1] || MemberQ[selectedSet, elem2],
        ret
    , With[{
        testRet = operator[elem1, elem2]
    }, If[testRet === None, ret, Append[ret, If[testRet =!= True, Append[pair, testRet], pair]]]]]]
], {}, Flatten[Table[{i, j}, {i, 1, Length@list}, {j, i + 1, Length@list}], 1]];
PdSymChristoffelToRiemann[expr_Plus, filter_] := Module[
    {exprL, vb, pdChrisList, chkList},
    exprL = MapIndexed[{#2[[1]], #1} &, List @@ expr];
    pdChrisList = Join @@ (ReplaceList[{
        HoldPattern[{n_, e: PD[-a0_?UpIndexQ]@chris_?SymChristoffelQ[d0_?UpIndexQ, -b0_?UpIndexQ, -c0_?UpIndexQ]}] :> {
            n,
            chris,
            {a0, b0, c0, d0},
            1
        } /; FilterExprList[filter, chris]
    ,
        HoldPattern[{n_, e: Times[
            PD[-a0_?UpIndexQ]@chris_?SymChristoffelQ[d0_?UpIndexQ, -b0_?UpIndexQ, -c0_?UpIndexQ], other__
        ]}] :> {
            n,
            chris,
            {a0, b0, c0, d0},
            Times[other]
        } /; FilterExprList[filter, chris]
    }] /@ exprL);
    chkList = FoldBinaryTest[TestPdSymChrisPair[exprL], pdChrisList];
    pdChrisList = With[{p1 = pdChrisList[[#1]], p2 = pdChrisList[[#2]]}, {p1[[1]], p2[[1]], PdSymChrisPairToSymChris2[p1, #3]}] & @@@ chkList;
    Plus @@ Join[Delete[exprL[[All, 2]], Transpose@{Flatten@pdChrisList[[All, 1 ;; 2]]}], pdChrisList[[All, 3]]]
];
PdSymChristoffelToRiemann[expr_, _] := expr;
PdSymChristoffelToRiemann[expr_] := PdSymChristoffelToRiemann[expr, All];
SyntaxInformation[PdSymChristoffelToRiemann] = {"ArgumentsPattern" -> {_, _.}};

ChangeCovDNonChristoffel[expr_, cd1_, cd2_] := expr /. HoldPattern[e: PD[_]@t_?NonChristoffelQ[___]] :> ChangeCovD[e, cd1, cd2];
SyntaxInformation[ChangeCovDNonChristoffel] = {"ArgumentsPattern" -> {_, _, _}};

TestCovDCommuPair[{_, cd_, {a1_, b1_}, term1_, other1_}, {_, cd_, {a2_, b2_}, term2_, other2_}] := If[ToCanonicalN[
    other1 * cd[b1]@cd[a1]@term1 + other2 * cd[a2]@cd[b2]@term2
] === 0, True, None];
CovDCommuPairToRiemann[{_, cd_, {a1_, b1_}, term1_, other1_}] := ToCanonicalN[-other1 * cd[b1]@cd[a1]@term1 + xAct`xTensor`Private`makeCommuteCovDs[term1, cd, {b1, a1}] * other1];
CovDCommuToRiemann[expr_Plus, filter_] := Module[
    {exprL, cdcdTerms},
    exprL = MapIndexed[{#2[[1]], #1} &, List @@ expr];
    cdcdTerms = Join @@ (ReplaceList[{
        HoldPattern[{n_, e: cd_?CovDQ[a1_]@cd_?CovDQ[a2_]@term_}] :> {n, cd, {a1, a2}, term, 1} /; FilterExprList[filter, cd],
        HoldPattern[{n_, e: Times[cd_?CovDQ[a1_]@cd_?CovDQ[a2_]@term_, other__]}] :> {n, cd, {a1, a2}, term, Times[other]} /; FilterExprList[filter, cd]
    }] /@ exprL);
    cdcdTerms = With[{p1 = cdcdTerms[[#1]], p2 = cdcdTerms[[#2]]}, {p1[[1]], p2[[1]], CovDCommuPairToRiemann@p1}] & @@@ FoldBinaryTest[TestCovDCommuPair, cdcdTerms];
    Plus @@ Join[Delete[exprL[[All, 2]], Transpose@{Flatten@cdcdTerms[[All, 1 ;; 2]]}], cdcdTerms[[All, 3]]]
];
CovDCommuToRiemann[expr_, _] := expr;
CovDCommuToRiemann[expr_] := CovDCommuToRiemann[expr, All];
SyntaxInformation[PdSymChristoffelToRiemann] = {"ArgumentsPattern" -> {_, _.}};

SignOfIndex[_?UpIndexQ] = 1;
SignOfIndex[-_?UpIndexQ] = -1;
SeparateMetricOneData[t_[inds__], ud_] := With[{
    li = MapThread[If[SignOfIndex@#1 =!= #2,
        With[{a1 = UniqueIndex@#1}, {-a1, delta[a1, #1]}]
    ,
        {#1, 1}
    ] &, {{inds}, ud}]
}, {t @@ li[[All, 1]], Times @@ li[[All, 2]]}];
SeparateMetricRiemann[expr_, filter_] := expr /. e: rt_?RiemannTensorQ[inds__] :> (Times @@ SeparateMetricOneData[e, {-1, -1, -1, -1}]) /; FilterExprList[filter, rt];
SeparateMetricRiemann[expr_] := SeparateMetricRiemann[expr, All];
SyntaxInformation[SeparateMetricRiemann] = {"ArgumentsPattern" -> {_, _.}};

ContractTensorWithMetric[expr_, filter_, metrics_] := expr //. {
    t_?xTensorQ[inds1__] metric_?MetricQ[inds2__] :> ContractMetric[t[inds1]metric[inds2], metrics] /; Length@Intersection[{inds1}, -{inds2}] === 2 && FilterExprList[filter, t]
};
SyntaxInformation[ContractTensorWithMetric] = {"ArgumentsPattern" -> {_, _, _}};

ListToCanonical[list_List, group_] := Module[
    {canonList, repeated, perm},
    canonList = Sort[list];
    repeated = Transpose[Position[canonList, #, {1}]][[1]] & /@ Union[list];
    perm = CanonicalPerm[Images@InversePermutation@Ordering@list, Length@list, group, {}, RepeatedSet /@ repeated];
    If[perm === 0, {0, list}, {If[Head[perm] === Times, -1, 1], PermuteList[canonList, InversePerm@perm]}]
];
SyntaxInformation[ListToCanonical] = {"ArgumentsPattern" -> {_, _}};

End[];

Protect @@ Names["`*"];

EndPackage[];