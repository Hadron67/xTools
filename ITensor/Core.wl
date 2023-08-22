BeginPackage["xTools`ITensor`Core`"]; (* depends on xAct`xPerm` secretly *)

Scan[Unprotect@#; ClearAll@#; &, Names@{"`*", "`Private`*"}];

(* interface *)
VSpaceQ::usage = "VSpaceQ[tslot]";
DimOfVSpace::usage = "";

(* IObject *)
ITensorMapSlot::usage = "ITensorMapSlot[bool]";
ITensorToScalarMapSlot::usage = "";
ITensorToTensorMapSlot::usage= "ITensorToTensorMapSlot[homoId]";
IObjectQ::usage = "";
SlotsOfIObject::usage = "SlotsOfIObject[obj[...]]";
SymmetryGroupOfIObject::usage = "SymmetryGroupOfIObject[obj] gives the symmetry group of the tensor slots of obj.";
CommutivityOfIObject::usage = "CommutivityOfIObject gives the commutivity of the tensor map slots of obj.";

Commutive::usage = "";
Anticommutive::usage = "";

(* IObject util *)
TensorMapGroups::usage = "";
ExtractAbsIndices::usage = "";
ApplyAbsIndices::usage = "";
ExtractScalarFunctionArgs::usage = "";
ApplyFunctionArgs::usage = "";
ExtractTensorFunctionArgs::usage = "";
ApplyTensorFunctionArgs::usage = "";

(* TODO: usages *)
AbsIndexQ::usage = "";
UpDownAbsIndexQ::usage = "";

MetricOfTensorSlot::usage = "MetricOfTensorSlot[slot] gives the ";

ScalarExpr::usage = "ScalarExpr[expr] akins to Scalar[expr] in xAct.";
FindIndices::usage = "";
FindFreeIndices::usage = "";
FindDummyIndices::usage = "";
FindFreeAndDummyIndices::usage = "";
CheckIndices::usage = "";

ISort::usage = "";
ISortedIObject::usage = "";
ISortedOther::usage = "";
IndexSlot::usage = "";

ToCanonical::usage = "";

(* utilities *)
NoneArg::usage = "";
DispatchIObjectArgs::usage = "";
NoneAnd::usage = "";

Begin["`Private`"];

AssertEvaluation::noeval = "Expression not evaluating: `1`.";
AssertEvaluation[expr_] := With[{val = expr}, If[Hold@expr === Hold@val, Message[AssertEvaluation::noeval, HoldForm@expr]; val, val]];
SetAttributes[AssertEvaluation, HoldAll];

(* default interface *)
IObjectQ[Integer|Rational|Complex|String|Real] = False;
IObjectQ[_Symbol] = True;
IObjectQ[Derivative[__]] = True;
SyntaxInformation[IObjectQ] = {"ArgumentsPattern" -> {_}};

SlotsOfIObject[Plus, args___] := ConstantArray[ITensorToTensorMapSlot@1, Length@{args}];
SlotsOfIObject[Derivative[__][fn_], args___] := SlotsOfIObject[fn, args];
SlotsOfIObject[fn_, args___] := ConstantArray[ITensorToTensorMapSlot@None, Length@{args}];
SyntaxInformation[SlotsOfIObject] = {"ArgumentsPattern" -> {__}};

SymmetryGroupOfIObject[Derivative[__][fn_], args___] := SymmetryGroupOfIObject[fn, args];
SymmetryGroupOfIObject[fn_, ___] := xAct`xPerm`StrongGenSet[{}, xAct`xPerm`GenSet[]];
SyntaxInformation[SymmetryGroupOfIObject] = {"ArgumentsPattern" -> {__}};

CommutivityOfIObject[b_Symbol, ___] := Commutive /; MemberQ[Attributes@b, Orderless];
CommutivityOfIObject[_, ___] = None;
SyntaxInformation@CommutivityOfIObject = {"ArgumentsPattern" -> {__}};

VSpaceQ[None|ITensorToScalarMapSlot|_ITensorToTensorMapSlot] = False;
VSpaceQ[_] = True;
SyntaxInformation[VSpaceQ] = {"ArgumentsPattern" -> {_}};

DimOfVSpace[a_] := a;
SyntaxInformation@DimOfVSpace = {"ArgumentsPattern" -> {_}};

TensorMapGroups[slots_] := With[{
    groups = GroupBy[Cases[MapIndexed[{#1, #2[[1]]} &, slots], {_ITensorToTensorMapSlot, _}], First@First@# &]
},
    Join[Thread@{Lookup[groups, None, {}][[All, 2]]}, Values[Delete[groups, Key@None]][[All, All, 2]]]
];
TensorMapGroups[slots_, args_] := args[[#]] & /@ TensorMapGroups@slots;
SyntaxInformation@TensorMapGroups = {"ArgumentsPattern" -> {_, _.}};

ExtractAbsIndices[slots_, args_] := Cases[Thread@{args, slots}, {_, _?VSpaceQ}][[All, 1]];
ExtractScalarFunctionArgs[slots_, args_] := Cases[Thread@{args, slots}, {_, ITensorToScalarMapSlot}][[All, 1]];
ExtractTensorFunctionArgs[slots_, args_] := Cases[Thread@{args, slots}, {_, _ITensorToTensorMapSlot}][[All, 1]];
SyntaxInformation@ExtractAbsIndices = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation@ExtractScalarFunctionArgs = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation@ExtractTensorFunctionArgs = {"ArgumentsPattern" -> {_, _}};

Null@NoneArg;
DispatchIObjectArgsOne[pat_][{ret_, args_}, slot_] := If[MatchQ[slot, pat], {Append[ret, args[[1]]], args[[2 ;;]]}, {Append[ret, NoneArg], args}];
DispatchIObjectArgs[slots_, pat_, args_] := Fold[DispatchIObjectArgsOne[pat], {{}, args}, slots][[1]];
SyntaxInformation@DispatchIObjectArgs = {"ArgumentsPattern" -> {_, _, _}};

NoneAnd[NoneArg..., a_, ___] := a /; a =!= NoneArg;
SyntaxInformation@NoneAnd = {"ArgumentsPattern" -> {___}};

AbsIndexQ[_Symbol] = True;
SyntaxInformation[AbsIndexQ] = {"ArgumentsPattern" -> {_}};

UpDownAbsIndexQ[_?AbsIndexQ] = True;
UpDownAbsIndexQ[-_?AbsIndexQ] = True;
SyntaxInformation@UpDownAbsIndexQ = {"ArgumentsPattern" -> {_}};

Options[FindIndices] = {CheckIndices -> False};
FindIndicesPlusLike[expr_, check_] := Union @@ (FindIndicesProductLike[#, check] & /@ expr);
FindIndicesProductLike[expr_, check_] := Union @@ (FindIndices /@ expr);
FindIndices[fn_?IObjectQ[args___]] := With[{
    slots = SlotsOfIObject[fn, args]
}, Union[FindIndicesPlusLike[TensorMapGroups[slots, {args}], OptionValue@CheckIndices], ExtractAbsIndices[slots, {args}]]];
FindIndices[ISortedIObject[_, _, _, dummies_, frees_]] := Join[dummies, frees];
FindIndices[__] = {};
SetAttributes[FindIndices, HoldFirst];
SyntaxInformation[FindIndices] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

DropPairs[inds_List] := Complement[inds, -inds];
TakePairs[inds_List] := Intersection[inds, -inds];
TakePairAbsIndices[inds_List] := TakePairs@Select[inds, UpDownAbsIndexQ];
DropPairAbsIndices[inds_List] := Complement[inds, -Select[inds, UpDownAbsIndexQ]];

Options@FindFreeIndices = Options@FindIndices;
FindFreeIndices[expr_, opt: OptionsPattern[]] := DropPairAbsIndices@FindIndices[expr, opt];
SyntaxInformation[FindFreeIndices] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

Options@FindDummyIndices = Options@FindIndices;
FindDummyIndices[expr_, opt: OptionsPattern[]] := Select[TakePairAbsIndices@FindIndices[expr, opt], AbsIndexQ];
SyntaxInformation@FindDummyIndices = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

Options@FindFreeAndDummyIndices = Options@FindIndices;
FindFreeAndDummyIndices[expr_, opt: OptionsPattern[]] := With[{inds = FindIndices[expr, opt]}, {DropPairAbsIndices@inds, Select[TakePairAbsIndices@inds, AbsIndexQ]}];
SyntaxInformation@FindFreeAndDummyIndices = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

ISort::uncomm = "Unknown communitivity `1`.";
ISort[expr: obj_?IObjectQ[args___]] := With[{
    slots = SlotsOfIObject[obj, args]
}, With[{
    sortedArgs = SortIObjectArgs[CommutivityOfIObject[obj, args], slots, {args}],
    inds = FindFreeAndDummyIndices@expr,
    symGroupOrder = AssertEvaluation@xAct`xPerm`OrderOfGroup@SymmetryGroupOfIObject[obj, args]
}, ISortedIObject[obj, sortedArgs, ##] &[symGroupOrder, Length@inds[[2]], Length@inds[[1]]]]];
ISort[expr_] := ISortedOther@expr;
SyntaxInformation[ISort] = {"ArgumentsPattern" -> {___}};

SetAttributes[{ISortedIObject, ISortedOther}, HoldAll];
IObjectQ@ISortedIObject ^= False;
IObjectQ@ISortedOther ^= False;
SyntaxInformation@ISortedIObject = {"ArgumentsPattern" -> {___}};
SyntaxInformation@ISortedOther = {"ArgumentsPattern" -> {___}};

SortIObjectArgs[None, _, args_] := args;
SortIObjectArgs[comm_, slots_, args_] := MapThread[NoneAnd, {
    DispatchIObjectArgs[slots, _ITensorToTensorMapSlot, SortIObjectList[comm, ISort /@ ExtractTensorFunctionArgs[slots, args]]],
    DispatchIObjectArgs[slots, ITensorToScalarMapSlot, SortIObjectList[comm, ISort /@ ExtractScalarFunctionArgs[slots, args]]],
    args
}];
SortIObjectList[Commutive, list_] := Sort[list, IObjectOrder];
SortIObjectList[comm_, _] := Null /; (Message[ISort::uncomm, comm]; False);

RemoveAbsIndices[slots_, args_] := DeleteCases[Thread@{args, slots}, {_, _?VSpaceQ}][[1]];
IObjectOrder[ISortedOther[e1_], ISortedOther[e2_]] := Order[e1, e2];
IObjectOrder[_ISortedIObject, _ISortedOther] = 1;
IObjectOrder[_ISortedOther, _ISortedIObject] = -1;
IObjectOrder[ISortedIObject[fn1_, {args1___}, order1__], ISortedIObject[fn2_, {args2___}, order2__]] := With[{
    i = -Order[{order1, fn1}, {order2, fn2}]
}, If[i === 0, Order[RemoveAbsIndices[SlotsOfIObject[fn1, args1], {args1}], RemoveAbsIndices[SlotsOfIObject[fn2, args2], {args2}]], i]];

ToCanonicalIObject[e_ISortedOther, ___] := e;
ToCanonicalIObject[ISortedIObject[obj_, args_, ___], opt: OptionsPattern[]] := Module[
    {},
    0
];

End[];

Protect @@ Names["`*"];

EndPackage[];
