BeginPackage["xTools`Misc`"];

(Unprotect[#]; ClearAll[#];) & /@ Names@{$Context <> "*", $Context <> "Private`*"};

FVariation::usage = "FVariation[expr, order] computes the functional variation.";
ConstantFunctions::usage = "ConstantFunctions is an option of FVariation that specifies constant functions.";
ConstantxTensors::usage = "ConstantxTensors is an option of FVariation that specifies constant tensors from xAct.";

DefPrintAs::usage = "DefPrintAs[name, string] makes name printed as the given string.";
SparseRowReduce::usage = "";
ReducedRowNumber::usage = "";

Begin["`Private`"];

Options[FVariation] = {ConstantFunctions -> {}};
FVariation[expr_] := FVariation[expr, 1];
FVariation[expr_, 0, opt___] := expr;
FVariation[expr_Plus, order_, opt___] := FVariation[#, order, opt] & /@ expr;
FVariation[expr_List, order_, opt___] := FVariation[#, order, opt] & /@ expr;
FVariation[xTools`xTension`ETensor[expr_, inds_], order_, opt___] := xTools`xTension`ETensor[
    FVariation[expr, order, opt],
    inds
];
FVariation[xTools`xDecomp`GCTensor[arr_, basis_], order_, opt___] := xTools`xDecomp`GCTensor[
    Map[FVariation[#, order, opt] &, arr, {Length@basis}],
    basis
];
FVariation[xTools`xDecomp`GCTensor[arr_, basis_][inds__], order_, opt___] := xTools`xDecomp`GCTensor[
    Map[FVariation[#, order, opt] &, arr, {Length@basis}],
    basis
][inds];
FVariation[expr_Times, order_Integer /; order >= 1, opt___] := FVariation[
    Plus @@ MapIndexed[FVariation[#, 1, opt]*Delete[expr, #2[[1]]] &, List @@ expr],
    order - 1,
    opt
];
FVariation[Power[a_, b_], order_Integer /; order >= 1, opt___] := FVariation[
    b*a^(b - 1) FVariation[a, 1, opt] + a^b Log[a] FVariation[b, 1, opt],
    order - 1,
    opt
];
FVariation[FVariation[expr_, order1_, opt___], order2_, opt___] := FVariation[expr, order1 + order2, opt];
FVariation[fn_Symbol[___], order_Integer /; order >= 1, opt: OptionsPattern[]] := 0 /; MemberQ[OptionValue[ConstantFunctions], fn];
FVariation[Derivative[__][fn_Symbol][__], order_Integer /; order >= 1, opt: OptionsPattern[]] := 0 /; MemberQ[OptionValue[ConstantFunctions], fn];
FVariation[_Symbol, order_Integer /; order >= 1, ___] = 0;
FVariation[_?NumberQ, order_Integer /; order >= 1, ___] = 0;
SyntaxInformation[FVariation] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

DefPrintAs[sym_, str_] := sym /: MakeBoxes[sym, StandardForm] := InterpretationBox[
    StyleBox[RowBox@{str}, AutoSpacing -> False, ShowAutoStyles -> False],
    sym,
    Editable -> False
];
SyntaxInformation[DefPrintAs] = {"ArgumentsPattern" -> {_, _}};

(* SparseRowReduce *)
FirstNonZeroPosition[list_List] := FirstPosition[list, n_ /; n =!= 0, {Length@list + 1}, {1}, Heads -> False][[1]];
FirstNonZeroPosition[arr_SparseArray] := With[{
    pos = Sort[arr["ExplicitPositions"]]
}, If[Length@pos === 0, Length@arr + 1, pos[[1, 1]]]];

NormalizeRow[list_] := With[{
    factor = list[[FirstNonZeroPosition[list]]]
}, If[NumberQ@factor, list / factor, Simplify[list / factor]]];
RemoveZeroFromSparseArray[arr_] := SparseArray[DeleteCases[ArrayRules[arr], _ -> 0], Dimensions@arr]; (* looks like SparseArray sometimes fails to remove zero elements *)
ResimplifyRow[list_List] := list;
ResimplifyRow[list_SparseArray] := RemoveZeroFromSparseArray@list;
SubstractRowWithNormalizedRow[row_, normalizedRow_] := With[{
    coef = row[[FirstNonZeroPosition@normalizedRow]]
}, If[coef =!= 0, ResimplifyRow@Simplify[row - normalizedRow * coef], row]];

SparseRowReduce[arr_List, row_Integer] := With[{
    chosenRow = row + Ordering[FirstNonZeroPosition /@ arr[[row ;;]], 1][[1]] - 1
}, If[FirstNonZeroPosition@arr[[chosenRow]] > Length@arr[[1]],
    arr,
    With[{
        newArr = If[row === chosenRow,
            ReplacePart[arr, {row -> NormalizeRow@arr[[row]]}],
            ReplacePart[arr, {row -> NormalizeRow@arr[[chosenRow]], chosenRow -> arr[[row]]}]
        ]
    }, With[{
        rowData = newArr[[row]]
    }, MapIndexed[If[#2[[1]] === row, rowData, SubstractRowWithNormalizedRow[#1, rowData]] &, newArr]]]
]];
SparseRowReduce[arr_SparseArray, row_] := SparseRowReduce[Extract[arr, Thread@{Range@Length@arr}], row];
SparseRowReduce[arr_, l_List] := Fold[SparseRowReduce, arr, l];
SparseRowReduce[arr_] := SparseRowReduce[arr, Range@Length@arr];

ReducedRowNumber[arr_] := With[{
    elemPos = FirstNonZeroPosition /@ arr,
    rows = Length@arr,
    cols = Length@arr[[1]]
}, NestWhile[# + 1 &, 1, # <= rows && With[{
    pos = elemPos[[#]]
}, If[pos > cols,
    AllTrue[elemPos[[# + 1 ;;]], # === pos &],
    arr[[#, pos]] === 1 && AllTrue[Range[# - 1], arr[[#, pos]] === 0 &] && AllTrue[elemPos[[# + 1 ;;]], # > pos &]
]] &]];

End[];

Protect @@ Names[$Context <> "*"];

EndPackage[];