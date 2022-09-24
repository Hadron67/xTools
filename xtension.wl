BeginPackage["xTools`xTension`", {"xAct`xTensor`"}];

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

DecomposedChristoffelRules::usage = "DecomposedChristoffelRules[cd] returns the decomposed Christoffel tensor rules of CovD cd, using cached result if already calculated.";

DecomposedRiemannRules::usage = "DecomposedRiemannRules[cd] returns the decomposed Riemann tensor rules of CovD cd, using cached result if already calculated.";

DecomposedRicciRules::usage = "DecomposedRicciRules[cd] returns the decomposed Ricci tensor rules of CovD cd, using cached result if already calculated.";

DecomposedRicciScalarRule::usage = "DecomposedRicciScalarRule[cd] returns the decomposed Ricci scalar rules of CovD cd, using cached result if already calculated.";

AllDecomposedRules::usage = "AllDecomposedRules[cd] returns all the decomposed rules of CovD cd, using cached result if already calculated.";

ClearDecomposedCache::usage = "ClearDecomposedCache[cd] clears all decomposed cache of Covd cd.";

ToDecomposed::usage = "ToDecomposed[expr, cd] or ToDecomposed[cd][expr] convert expr into decomposed forms of lists using rules AllDecomposedRules[cd].";

ReplaceIndicesRules::usage = "ReplaceIndicesRules[expr, fromVB, toVB, opt] returns a list of rules that replaces all indices of fromVB in expr to toVB.";

ReplaceCovDs::usage = "ReplaceCovDs[expr, cdFrom, cdTo, opt] replaces all indices and associated tensors (curvatures, etc) of cdFrom in expr with that of cdTo.";

OtherRules::usage = "OtherRules is an option for ReplaceIndicesRules and ReplaceCovDs that defines other rules to be applied.";

Protect[OtherRules];

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
Protect@MakeDecomposedRules;

SetDecomposedRules[tensor_[inds__], vals_] := (#1 := #2) & @@@ MakeDecomposedRules[tensor[inds],vals];
Protect@SetDecomposedRules;

BasisOfCoordinateParameter[_]:=Null;
DualBasisOfCoordinateParameter[_]:=Null;
CoordinateParameterOfVBundle[_]:=Null;
VBundleOfCoordinateParameter[_]:=Null;
CoordinateBasisQ[_]:=False;
Protect[BasisOfCoordinateParameter, DualBasisOfCoordinateParameter, CoordinateParameterOfVBundle, VBundleOfCoordinateParameter, CoordinateBasisQ];

DefCoordinateParameter::nsvb="`1` is not a subvbundle of `2`";
DefCoordinateParameter[parentVB_?VBundleQ -> vb_?VBundleQ, param_?ParameterQ, e_, ed_]:=Module[
    {a, b, c, otherVB},
    If[parentVB == vb || !SubvbundleQ[parentVB, vb], Throw[Message[DefCoordinateParameter::nsvb, vb, parentVB]]];

    {a} = GetIndicesOfVBundle[vb, 1];
    otherVB = DeleteCases[List @@ First@SplittingsOfVBundle@parentVB, vb];
    DefTensor[e[a], BaseOfVBundle@vb, PrintAs -> "(\!\(\*SubscriptBox[\(\[PartialD]\), \(" <> PrintAs[param] <> "\)]\))"];
    DefTensor[ed[-a], BaseOfVBundle@vb, PrintAs->"(d" <> PrintAs[param] <> ")"];
    e /: e[b_Symbol] ed[-b_Symbol] := 1;
    CoordinateBasisQ[e] ^= True;
    CoordinateBasisQ[ed] ^= True;

    With[
        {pmQ = Symbol[ToString[#] <> "`pmQ"]},
        param /: PD[_?pmQ]@param:=0;
        e /: PD[_?pmQ]@e[_Symbol]:=0;
        ed /: PD[_?pmQ]@ed[-_Symbol]:=0;
    ] & /@ otherVB;
    With[
        {cd = CovDOfMetric@#, pmQ = Symbol[ToString[VBundleOfMetric@#] <> "`pmQ"]},
        param /: cd[_?pmQ]@param := 0;
        e /: cd[_?pmQ]@e[_Symbol] := 0;
        ed /: cd[_?pmQ]@ed[-_Symbol] := 0;
    ] & /@ Select[Union @@ MetricsOfVBundle /@ otherVB, CovDOfMetric[#] =!= PD &];

    CoordinateParameterOfVBundle[vb] ^= param;
    VBundleOfCoordinateParameter[param] ^= vb;
    BasisOfCoordinateParameter[param] ^= e;
    DualBasisOfCoordinateParameter[param] ^= ed;
];
Protect@DefCoordinateParameter;

DropCoordinateBasis[expr_] := expr /. _?CoordinateBasisQ[_] -> 1;
CoordinateParameterQ[p_] := BasisOfCoordinateParameter[p] =!= Null;
CoordinateParameterOfIndex[i_] := CoordinateParameterOfVBundle@VBundleOfIndex@UpIndex@i;
CoordinateIndexQ[i_] := AIndexQ[i] && CoordinateParameterOfIndex@i =!= Null;
Protect[DropCoordinateBasis, CoordinateParameterQ, CoordinateParameterOfIndex, CoordinateIndexQ];

Unprotect@{PD,ParamD};
PD[-a_Symbol?CoordinateIndexQ][A_] := With[
    {param = CoordinateParameterOfIndex[a]},
    ParamD[param][A] DualBasisOfCoordinateParameter[param][-a]
];
ParamD[p_Symbol?CoordinateParameterQ]@PD[-a_Symbol]@A_ := PD[-a]@ParamD[p]@A;
Protect@{PD,ParamD};

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
SplitAllIndices[expr_, vb_] := expr // SplitIndex[#, DecomposedIndicesList[FindFreeIndices@Evaluate@#, vb]] & // ScalarTraceProductDummy // TraceProductDummy;
SplitAllIndices[vb_][expr_] := SplitAllIndices[expr, vb];
Protect@SplitAllIndices;

ToCanonicalN[expr_, opt___] := ToCanonical[expr, UseMetricOnVBundle -> None, opt];
SimplificationN[expr_] := Implode[expr, ParamD] // Simplification // Explode;
Protect[ToCanonicalN, SimplificationN];

CalculateDecomposed[expr_, vb_, calc_, calc2_] := Module[
    {freeInds, indsMap, vals, lhs},
    freeInds = FindFreeIndices@Evaluate@expr;
    indsMap = DecomposedIndicesList[freeInds, vb];
    vals = expr // calc // SplitIndex[#, indsMap] & // TraceProductDummy // Flatten // calc2 // SimplificationN;
    lhs = expr // SplitIndex[#, indsMap] & // Flatten;
    Flatten[MakeRule[Evaluate@#, MetricOn -> None] & /@ Transpose@{lhs, vals}]
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
Protect[SelectNonFlatCD, CalculateDecomposedChristoffel];

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
Protect@CalculateDecomposedRiemann;

CalculateDecomposedRicci[cd_, riemannRules_] := Module[
    {metric, riemann, ricci, a, b, c, d},
    metric = MetricOfCovD@cd;
    riemann = Riemann@cd;
    ricci = Ricci@cd;
    {a, b, c, d} = GetIndicesOfVBundle[VBundleOfMetric@metric, 4];
    CalculateDecomposed[ricci[-a, -b], VBundleOfMetric@metric, riemann[-a, -c, -b, -d] metric[c, d] &, e |-> e /. riemannRules]
];
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
    MakeRule[{Evaluate[RicciScalar[cd][]], Evaluate@Scalar@ret}, 
    MetricOn -> None]
];
Protect@CalculateDecomposedRicciScalar;

DecomposedChristoffelRulesOf[_] := Null;
DecomposedRiemannRulesOf[_] := Null;
DecomposedRicciRulesOf[_] := Null;
DecomposedRicciScalarRuleOf[_] := Null;
Protect[DecomposedChristoffelRulesOf, DecomposedRiemannRulesOf, DecomposedRicciRulesOf, DecomposedRicciScalarRuleOf];

Lazy[expr_, cache_] := If[cache === Null, cache ^= expr, cache];
SetAttributes[Lazy, HoldAll];
Protect@Lazy;

DecomposedChristoffelRules[cd_] := Lazy[CalculateDecomposedChristoffel[cd], DecomposedChristoffelRulesOf[cd]];
DecomposedRiemannRules[cd_] := Lazy[CalculateDecomposedRiemann[cd, DecomposedChristoffelRules[cd]], DecomposedRiemannRulesOf[cd]];
DecomposedRicciRules[cd_] := Lazy[CalculateDecomposedRicci[cd, DecomposedRiemannRules[cd]], DecomposedRicciRulesOf[cd]];
DecomposedRicciScalarRule[cd_] := Lazy[CalculateDecomposedRicciScalar[cd, DecomposedRicciRules[cd]], DecomposedRicciScalarRuleOf[cd]];
Protect[DecomposedChristoffelRules, DecomposedRiemannRules, DecomposedRicciRules, DecomposedRicciScalarRule];
AllDecomposedRules[cd_] := Flatten@{
    DecomposedChristoffelRules[cd],
    DecomposedRiemannRules[cd],
    DecomposedRicciRules[cd],
    DecomposedRicciScalarRule[cd]
};
ClearDecomposedCache[cd_] := (
    cd /: DecomposedChristoffelRulesOf[cd] =.;
    cd /: DecomposedRiemannRulesOf[cd] =.;
    cd /: DecomposedRicciRulesOf[cd] =.;
    cd /: DecomposedRicciScalarRuleOf[cd] =.;
);
Protect[AllDecomposedRules, ClearDecomposedCache];

ScalarTraceProductDummy[expr_] := expr /. Scalar[e_] :> Scalar@TraceProductDummy@e;
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
        // SplitIndex[#, DecomposedIndicesList[FindFreeIndices@Evaluate@#, VBundleOfMetric@metric]] &
        // ScalarTraceProductDummy
        // TraceProductDummy
        // ReplaceAll@AllDecomposedRules[cd]
        // Fold[ChangeCovD[#1, PD, #2] &, #, cds] &
];
Protect[ScalarTraceProductDummy, ToDecomposed];

Options[ReplaceIndicesRules] = {OtherRules -> {}};
ReplaceIndicesRules::nsi = "index(ices) `1` not found in expression";
ReplaceIndicesRules[expr_, fromVB_, toVB_, opt: OptionsPattern[]] := Module[
    {inds, extraRules, notFoundInds},
    extraRules = OptionValue[OtherRules];
    inds = UpIndex /@ Select[List @@ FindIndices@expr, VBundleOfIndex[#] === fromVB &] // DeleteDuplicates;
    notFoundInds = DeleteCases[#1 & @@@ extraRules, Alternatives @@ inds];
    If[Length@notFoundInds > 0, Message[ReplaceIndicesRules::nsi, notFoundInds]];
    inds = DeleteCases[inds, Alternatives @@ (#1 & @@@ extraRules)];
    Join[extraRules, #1 -> #2 & @@@ Transpose@{inds, GetIndicesOfVBundle[toVB, Length@inds, #2 & @@@ extraRules]}]
];
Protect@ReplaceIndicesRules;

Options[ReplaceCovDs] = {OtherRules -> {}};
ReplaceCovDs[expr_, cdFrom_, cdTo_, opt: OptionsPattern[]] := 
  expr /. Flatten@{
    Select[OptionValue[OtherRules], ! AIndexQ[#[[1]]] &],
    ReplaceIndicesRules[
        expr,
        First@VBundlesOfCovD@cdFrom,
        First@VBundlesOfCovD@cdTo,
        OtherRules -> Select[OptionValue[OtherRules], AIndexQ[#[[1]]] &]
    ],
    cdFrom -> cdTo,
    MetricOfCovD[cdFrom] -> MetricOfCovD[cdTo],
    Riemann[cdFrom] -> Riemann[cdTo],
    Ricci[cdFrom] -> Ricci[cdTo],
    RicciScalar[cdFrom] -> RicciScalar[cdTo]
};
Protect@ReplaceCovDs;

End[];

EndPackage[];