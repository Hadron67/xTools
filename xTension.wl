BeginPackage["xTools`xTension`", {"xAct`xCore`", "xAct`xTensor`"}];

(Unprotect[#]; Remove[#];) & /@ Names@{"xTools`xTension`*", "xTools`xTension`Private`*"};

ToCanonicalN::usage = "Calls ToCanonical with UseMetricOnVBundle -> None.";
SimplificationN::usage = "Calls Simplification[] using Implode.";

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

SetCMetricRule::usage = "SetCMetricRule[metric, cmetric] defines component rules for metric using CTensor metric cmetric.";

UnsetCMetricRule::usage = "UnsetCMetricRule[metric] removes all rules defined by SetCMetricRule[metric, ...].";

WithCMetricRule::usage = "WithCMetricRule[expr_, metric_, cmetric_] evaluates expr by calling SetCMetricRule[metric, cmetric] first. It's also an option for ReplaceCovDs determining whether wrap the whole expression with WithCMetricRule or not.";

BHTemperature::usage = "BHTemperature[gtt, grr, gxx, r, r0] calculates the Hawking temperature of black hole with radial coordinate r and metric components gtt, grr, gxx, and horizon r0.";

OtherRules::usage = "OtherRules is an option for ReplaceIndicesRules and ReplaceCovDs that defines other rules to be applied.";

MetricInv::usage = "MetricInv is an option for SetCMetricRule that specifies the inverse metric explicitly.";

DefMetricNsd::usage = "DefMetricNsd[metric[-a, -b], covd, ...] calls DefMetric[1, metric[-a, -b], covd, ...] and then defines SignDetOfMetric[metric] as poison.";

FastReplaceDummies::usage = "FastReplaceDummies[expr] is a faster implementation of ReplaceDummies[expr].";

ETensor::usage = "ETensor[expr, {freeIndices}] represents an expression with dummy and free indices.";

ETensorProduct::usage = "ETensorProduct[e1, e2, ...] computes the tensor outer product of ETensors e1, e2, ...";

ETensorTranspose::usage = "ETensorTranspose[ETensor[...], perms] performs the transpose of the ETensor.";

ETensorDot::usage = "ETensorDot[T1...] computes the tensor dot for tensors expressed in ETensor form.";

ETensorPD::usage = "ETensorPD[T, vb] or ETensorPD[vb][T] adds a PD to the ETensor.";

ETensorContractTwo::usage = "ETensorContractTwo[T1, T2, {n1...}, {n2...}] contracts the n1... axes of T1 to n2.. axes of T2.";

ETensorRank::usage = "ETensorRank[ETensor[...]] gives the rank of the ETensor.";

UniqueIndex::usage = "UniqueIndex[a] gets a unique and temporary a-index of the form a$*.";

Protect[OtherRules, MetricInv];

Begin["`Private`"];

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
Protect@MakeDecomposedRules;

SetDecomposedRules[tensor_[inds__], vals_] := (#1 := #2) & @@@ MakeDecomposedRules[tensor[inds], vals];
SyntaxInformation[SetDecomposedRules] = {"ArgumentsPattern" -> {_, _}};
Protect@SetDecomposedRules;

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
Protect[BasisOfCoordinateParameter, DualBasisOfCoordinateParameter, CoordinateParameterOfVBundle, VBundleOfCoordinateParameter, CoordinateBasisQ];

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
Protect[DefCoordinateParameter, DropCoordinateBasis];

CoordinateParameterQ[p_] := BasisOfCoordinateParameter[p] =!= Null;
CoordinateParameterOfIndex[i_] := CoordinateParameterOfVBundle@VBundleOfIndex@UpIndex@i;
CoordinateIndexQ[i_] := AIndexQ[i] && CoordinateParameterOfIndex@i =!= Null;
Protect[CoordinateParameterQ, CoordinateParameterOfIndex, CoordinateIndexQ];

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
Protect@DecomposedIndicesList;
ScalarTraceProductDummy[expr_] := expr /. Scalar[e_] :> Scalar@TraceProductDummy@e;
Protect@ScalarTraceProductDummy;
SplitAllIndices[l_List, vb_] := SplitAllIndices[#, vb] & /@ l;
SplitAllIndices[expr_, vb_] := expr // SplitIndex[#, DecomposedIndicesList[FindFreeIndices@#, vb]] & // ScalarTraceProductDummy // TraceProductDummy;
SplitAllIndices[vb_][expr_] := SplitAllIndices[expr, vb];
SyntaxInformation[SplitAllIndices] = {"ArgumentsPattern" -> {_, _.}};
Protect@SplitAllIndices;

ToCanonicalN[expr_, opt___] := ToCanonical[expr, UseMetricOnVBundle -> None, opt];
SimplificationN[expr_] := Implode[expr, ParamD] // Simplification // Explode;
Protect[ToCanonicalN, SimplificationN];

CalculateDecomposed[expr_, vb_, calc_, calc2_] := Module[
    {freeInds, indsMap, vals, lhs},
    freeInds = FindFreeIndices@expr;
    indsMap = DecomposedIndicesList[freeInds, vb];
    vals = expr // calc // SplitIndex[#, indsMap] & // TraceProductDummy // Flatten // calc2 // SimplificationN;
    lhs = expr // SplitIndex[#, indsMap] & // Flatten;
    Flatten[MakeRule[#, MetricOn -> None] & /@ Transpose@{lhs, vals}]
];
Protect@CalculateDecomposed;

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
Protect[SelectNonFlatCD, CalculateDecomposedChristoffel];

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
Protect@ChristoffelToRiemann;

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
Protect@CalculateDecomposedRiemann;

CalculateDecomposedRicci[cd_, riemannRules_] := Module[
    {metric, riemann, ricci, a, b, c, d},
    metric = MetricOfCovD@cd;
    riemann = Riemann@cd;
    ricci = Ricci@cd;
    {a, b, c, d} = GetIndicesOfVBundle[VBundleOfMetric@metric, 4];
    CalculateDecomposed[ricci[-a, -b], VBundleOfMetric@metric, riemann[-a, -c, -b, -d] metric[c, d] &, e |-> e /. riemannRules]
];
SyntaxInformation[CalculateDecomposedRicci] = {"ArgumentsPattern" -> {_, _}};
Protect@CalculateDecomposedRicci;

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
Protect@CalculateDecomposedRicciScalar;

DecomposedChristoffelRules[cd_] := DecomposedChristoffelRules[cd] ^= CalculateDecomposedChristoffel[cd];
DecomposedRiemannRules[cd_] := DecomposedRiemannRules[cd] ^= CalculateDecomposedRiemann[cd, DecomposedChristoffelRules[cd]];
DecomposedRicciRules[cd_] := DecomposedRicciRules[cd] ^= CalculateDecomposedRicci[cd, DecomposedRiemannRules[cd]];
DecomposedRicciScalarRule[cd_] := DecomposedRicciScalarRule[cd] ^= CalculateDecomposedRicciScalar[cd, DecomposedRicciRules[cd]];
SyntaxInformation[DecomposedChristoffelRules] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[DecomposedRiemannRules] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[DecomposedRicciRules] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[DecomposedRicciScalarRule] = {"ArgumentsPattern" -> {_}};
Protect[DecomposedChristoffelRules, DecomposedRiemannRules, DecomposedRicciRules, DecomposedRicciScalarRule];
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
Protect[AllDecomposedRules, ClearDecomposedCache];

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
Protect@ToDecomposed;

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
Protect@ReplaceIndicesRules;

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
Protect@ReplaceCovDs;

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
Protect[DefRiemannVarD, UndefRiemannVarD];

Options[SetCMetricRule] = {MetricInv -> None};
SetCMetricRule[metric_, cmetric_, opt: OptionsPattern[]] := With[{
    inv = If[# =!= None, #, Inv@cmetric] &@OptionValue@MetricInv
},
    metric[a_?DownIndexQ, b_?DownIndexQ] := cmetric[a, b];
    metric[a_?UpIndexQ, b_?UpIndexQ] := inv[a, b];
];
SyntaxInformation[SetCMetricRule] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
Protect@SetCMetricRule;

UnsetCMetricRule[metric_] := (
    metric[a_?DownIndexQ, b_?DownIndexQ] =.;
    metric[a_?UpIndexQ, b_?UpIndexQ] =.;
);
SyntaxInformation[UnsetCMetricRule] = {"ArgumentsPattern" -> {_}};
Protect@UnsetCMetricRule;

Options[WithCMetricRule] = Options[SetCMetricRule];
WithCMetricRule[expr_, metric_, cmetric_, opt: OptionsPattern[]] := WithCleanup[
    SetCMetricRule[metric, cmetric, opt],
    expr,
    UnsetCMetricRule[metric]
];
SetAttributes[WithCMetricRule, HoldFirst];
SyntaxInformation[WithCMetricRule] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};
Protect@WithCMetricRule;

BHTemperature[gtt_, grr_, gxx_, r_, r0_] := With[{fp = D[gtt/gxx, r], hp = D[1/grr, r]}, 1/(4 Pi) Sqrt[gxx fp hp] // Simplify // ReplaceAll[r -> r0]];
SyntaxInformation[BHTemperature] = {"ArgumentsPattern" -> {_, _, _, _, _}};
Protect@BHTemperature;

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
Protect[DefMetricNsd];

FastReplaceDummies[expr_] := With[{
    rep = # -> UniqueIndex@# & /@ Union[List @@ UpIndex /@ FindDummyIndices[expr]]
}, expr /. rep];
SyntaxInformation[FastReplaceDummies] = {"ArgumentsPattern" -> {_}};
Protect[FastReplaceDummies];

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
ETensor[expr_] := ETensor[expr, List @@ FindFreeIndices@expr];
ETensor[expr_, {}] := Scalar[expr];
ETensor[expr_, inds_List][inds2__] /; Length@inds === Length@{inds2} := Module[
    {indPairs, matchedInds, unmatchedInds, unmatchedDummies, deltas},
    indPairs = Thread@{inds, {inds2}};
    matchedInds = Select[indPairs, SameIndexUDQ[#[[1]], #[[2]]] &];
    unmatchedInds = Select[indPairs, !SameIndexUDQ[#[[1]], #[[2]]] &];
    unmatchedDummies = UniqueIndex[UpIndex[#[[1]]]] & /@ unmatchedInds;
    unmatchedInds = Append[#1, #2] & @@@ Transpose@{unmatchedInds, unmatchedDummies};
    deltas = Times @@ (delta[ChangeIndex[CopyIndexSign[#1, #3]], #2] & @@@ unmatchedInds);
    (FastReplaceDummies[expr] /. (UpIndex[#1] -> #3 & @@@ unmatchedInds))*deltas /. (#1 -> #2 & @@@ matchedInds)
];
ETensor /:
    ETensor[expr1_, inds1_List] + ETensor[expr2_, inds2_List] /; CompatibleIndexListsQ[inds1, inds2] := (
        ETensor[
            FastReplaceDummies[expr1]
            + (FastReplaceDummies[expr2] /. Thread[inds2 -> inds1])
        , inds1]
    );
ETensor /: Times[ETensor[expr_, inds_], factors__] := (
    If[Cases[{factors}, _ETensor], Message[ETensor::invldmlt]];
    ETensor[Times[expr, factors], inds]
);
ETensor /: ToCanonical[ETensor[expr_, args__], opt___] := ETensor[ToCanonical[expr, opt], args];
ETensor /: Simplification[ETensor[expr_, args__], opt___] := ETensor[Simplification[expr, opt], args];
ETensor /: Simplify[ETensor[expr_, args__], opt___] := ETensor[Simplify[expr, opt], args];
ETensor /: ParamD[params__][ETensor[expr_, args__]] := ETensor[ParamD[params][expr], args];
ETensor /: FindFreeIndices[ETensor[_, inds_]] := IndexList @@ inds;
ETensor /: ScreenDollarIndices[ETensor[expr_, inds_]] := Module[
    {dollars, rep},
    dollars = Select[Union[UpIndex /@ inds], MemberQ[Attributes@#, Temporary] &];
    rep = If[Length@dollars > 0, Thread[dollars -> GetIndicesOfVBundle[VBundleOfIndex@UpIndex@First@dollars, Length@dollars, UpIndex /@ FindIndices@expr]], {}];
    ETensor[ScreenDollarIndices[expr /. rep], inds /. rep]
];
ETensor /: ChangeCovD[ETensor[expr_, inds_], args___] := ETensor[ChangeCovD[expr, args], inds];

SyntaxInformation[ETensor] = {"ArgumentsPattern" -> {_, _.}};
Protect[ETensor];

DedupeRules[inds1_, inds2_] := With[{
    dupes = Intersection[UpIndex /@ inds1, UpIndex /@ inds2]
}, If[Length@dupes > 0, Thread[dupes -> (UniqueIndex /@ dupes)], {}]];

IndexedScalarQ[_ETensor | _List] = False;
IndexedScalarQ[_] = True;

ETensorProduct[ETensor[expr1_, inds1_], ETensor[expr2_, inds2_], rest___] := With[{
    rep = DedupeRules[inds1, inds2]
},
    ETensorProduct[
        ETensor[FastReplaceDummies[expr1] * (FastReplaceDummies[expr2] /. rep), Join[inds1, inds2 /. rep]],
        rest
    ]
];
ETensorProduct[left___, Zero, right___] = Zero;
ETensorProduct[e_] := e;
ETensorProduct[expr_ /; !MatchQ[expr, _List], l_List, rest___] := ETensorProduct[expr, #, rest] & /@ l;
ETensorProduct[l_List, rest___] := ETensorProduct[#, rest] & /@ l;
ETensorProduct[ETensor[expr_, inds_], x_?IndexedScalarQ, rest___] := ETensorProduct[ETensor[expr * x, inds], rest];
ETensorProduct[x_?IndexedScalarQ, ETensor[expr_, inds_], rest___] := ETensorProduct[ETensor[expr * x, inds], rest];
ETensorProduct[x_?IndexedScalarQ, y_?IndexedScalarQ, rest___] := x*y*If[Length@{rest} === 0, 1, ETensorProduct[rest]];
SyntaxInformation[ETensorProduct] = {"ArgumentsPattern" -> {___}};
Protect[ETensorProduct];

ETensorDot0[ETensor[expr1_, {left___, a_}], ETensor[expr2_, {b_, right___}]] := With[{
    bi = UniqueIndex@b
},
    ETensor[FastReplaceDummies[expr1] * (FastReplaceDummies[expr2] /. b -> bi) * delta[ChangeIndex@a, ChangeIndex@bi], {left, right}]
];

ETensorDot[a_ETensor, b_ETensor, rest___] := ETensorDot[ETensorDot0[a, b /. DedupeRules[a[[2]], b[[2]]]], rest];

ETensorDot[left___, Zero, right___] = Zero;
ETensorDot[ETensor[expr_, inds_], x_?IndexedScalarQ, rest___] := ETensorDot[ETensor[expr * x, inds], rest];
ETensorDot[x_?IndexedScalarQ, ETensor[expr_, inds_], rest___] := ETensorDot[ETensor[expr * x, inds], rest];

ETensorDot[expr_ /; !MatchQ[expr, _List], l_List, rest___] := ETensorDot[expr, #, rest] & /@ l;
ETensorDot[l_List, rest___] := ETensorDot[#, rest] & /@ l;

ETensorDot[e_] := e;
SyntaxInformation[ETensorDot] = {"ArgumentsPattern" -> {___}};
Protect[ETensorDot];

ETensorContractTwo0[expr1_, inds1_, expr2_, inds2_, n1_List, n2_List] := With[{
    a = inds1[[#]] & /@ n1,
    b = inds2[[#]] & /@ n2
},  With[{bi = UniqueIndex /@ b},
    ETensor[FastReplaceDummies[expr1] * (FastReplaceDummies[expr2] /. Thread[b -> bi]) * (Times @@ (delta[ChangeIndex@#1, ChangeIndex@#2] & @@@ Thread@{a, bi})), Join[Delete[inds1, Transpose@{n1}], Delete[inds2, Transpose@{n2}]]]
]];

ETensorContractTwo[ETensor[expr1_, inds1_], ETensor[expr2_, inds2_], n1_, n2_] := With[{
    rep = DedupeRules[inds1, inds2]
}, ETensorContractTwo0[expr1, inds1, expr2 /. rep, inds2 /. rep, n1, n2]];
ETensorContractTwo[ETensor[expr_, inds_], x_?IndexedScalarQ, {}, {}] := ETensor[expr * x, inds];
ETensorContractTwo[x_?IndexedScalarQ, ETensor[expr_, inds_], {}, {}] := ETensor[expr * x, inds];
ETensorContractTwo[x_?IndexedScalarQ, y_?IndexedScalarQ, {}, {}] := x * y;
SyntaxInformation[ETensorContractTwo] = {"ArgumentsPattern" -> {_, _, _, _}};
Protect[ETensorContractTwo];

ETensorRank[ETensor[_, inds_]] := Length@inds;
ETensorRank[_] = 0;
SyntaxInformation[ETensorRank] = {"ArgumentsPattern" -> {_}};
Protect[ETensorRank];

ETensorTranspose[ETensor[expr_, inds_], perms_] := ETensor[expr, Permute[inds, perms]];
ETensorTranspose[expr_, {}] := expr;
ETensorTranspose[expr_, {_Integer}] := expr;
ETensorTranspose[0, _] = 0;
ETensorTranspose[Zero, _] = Zero;
ETensorTranspose[expr_Plus, perms_] := ETensorTranspose[#, perms] /@ expr;
SyntaxInformation[ETensorTranspose] = {"ArgumentsPattern" -> {_, _}};
Protect[ETensorTranspose];

ETensorPD[vb_][expr_] := ETensorPD[expr, vb];
ETensorPD[x_?IndexedScalarQ, vb_] := With[{
    ai = UniqueIndex@First@GetIndicesOfVBundle[vb, 1]
}, ETensor[PD[-ai]@x, {-ai}]];
ETensorPD[ETensor[expr_, inds_], vb_] := With[{
    ai = UniqueIndex@First@GetIndicesOfVBundle[vb, 1]
}, ETensor[PD[-ai]@expr, Prepend[inds, -ai]]];
SyntaxInformation[ETensorPD] = {"ArgumentsPattern" -> {_, _}};
Protect[ETensorPD];

End[];

EndPackage[];