BeginPackage["xTools`Misc`"];

(Unprotect[#]; ClearAll[#];) & /@ Names@{$Context <> "*", $Context <> "Private`*"};

FVariation::usage = "FVariation[expr, order] computes the functional variation.";
ConstantFunctions::usage = "ConstantFunctions is an option of FVariation that specifies constant functions.";
ConstantxTensors::usage = "ConstantxTensors is an option of FVariation that specifies constant tensors from xAct.";

DefPrintAs::usage = "DefPrintAs[name, string] makes name printed as the given string.";
SparseRowReduce::usage = "";
RowMultiplier;
RowNormalizer;
RowSimplifier;
MaxColumn::usage = "";
ReducedRowNumber::usage = "";
DropTraillingZeroRows::usage = "";

TermToList::usage = "TermToList[term]";
SeparateFactor::usage = "SeparateFactor[expr, fns]";
GroupPolynomialBy::usage = "GroupPolynomialBy[expr, fn]";
PolynomialToVec::usage = "PolynomialToVec[expr, fn, terms]";
CollectBy::usage = "CollectBy[expr, fn, action]";

SaveNotebookData::usage = "SaveNotebookData[names]";
LoadNotebookData::usage = "LoadNotebookData[]";
SavableObjQ::usage = "SavableObjQ[symbol]";

WriteDefinitionsAsDelayedValueTo::usage = "SaveAsSetDelayed[stream, symbol]";

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

NormalizeRow[list_, multiplier_, normalizer_] := With[{
    factor = normalizer@list[[FirstNonZeroPosition[list]]]
}, If[NumberQ@factor, multiplier[factor, list], Simplify[multiplier[factor, list]]]];
RemoveZeroFromSparseArray[arr_] := SparseArray[DeleteCases[ArrayRules[arr], _ -> 0], Dimensions@arr]; (* looks like SparseArray sometimes fails to remove zero elements *)
ResimplifyRow[list_List] := list;
ResimplifyRow[list_SparseArray] := RemoveZeroFromSparseArray@list;
SubstractRowWithNormalizedRow[row_, normalizedRow_, multiplier_, simplifier_] := With[{
    coef = row[[FirstNonZeroPosition@normalizedRow]]
}, If[coef =!= 0, ResimplifyRow@simplifier[row - multiplier[coef, normalizedRow]], row]];

Options[SparseRowReduce] = {
    MaxColumn -> Max,
    ProgressReporting -> Nothing,
    RowMultiplier -> Times,
    RowNormalizer -> (1/# &),
    RowSimplifier -> Simplify
};
SparseRowReduce[arr_List, row_Integer, opt : OptionsPattern[]] := With[{
    chosenRow = row + Ordering[FirstNonZeroPosition /@ arr[[row ;;]], 1][[1]] - 1,
    maxColumn = With[{op = OptionValue@MaxColumn}, If[op === Max, Length@arr[[1]], If[op < 0, op + Length@arr[[1]], op]]],
    rowMultiplier = OptionValue@RowMultiplier,
    rowNormalizer = OptionValue@RowNormalizer,
    rowSimplifier = OptionValue@RowSimplifier
}, WithCleanup[
    If[FirstNonZeroPosition@arr[[chosenRow]] > maxColumn,
        arr,
        With[{
            newArr = If[row === chosenRow,
                ReplacePart[arr, {row -> NormalizeRow[arr[[row]], rowMultiplier, rowNormalizer]}],
                ReplacePart[arr, {row -> NormalizeRow[arr[[chosenRow, rowMultiplier, rowNormalizer]]], chosenRow -> arr[[row]]}]
            ]
        }, With[{
            rowData = newArr[[row]]
        }, MapIndexed[If[#2[[1]] === row, rowData, SubstractRowWithNormalizedRow[#1, rowData, rowMultiplier, rowSimplifier]] &, newArr]]]
    ],
    OptionValue[ProgressReporting][row, Length@arr]
]];
SparseRowReduce[arr_SparseArray, row_, opt___] := SparseRowReduce[Extract[arr, Thread@{Range@Length@arr}], row, opt];
SparseRowReduce[arr_, l_List, opt___] := Fold[SparseRowReduce[#1, #2, opt] &, arr, l];
SparseRowReduce[arr_, All, opt___] := SparseRowReduce[arr, Range@Length@arr, opt];
SparseRowReduce[arr_] := SparseRowReduce[arr, Range@Length@arr];
SyntaxInformation@SparseRowReduce = {_, _, OptionsPattern[]};

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

AllZeroQ[list_List] := AllTrue[list, # === 0 &];
AllZeroQ[arr_SparseArray] := arr["ExplicitPositions"] === {};

DropTraillingZeroRows[list_List] := Drop[list, -Length@list + NestWhile[# - 1 &, Length@list, # >= 0 && AllZeroQ[list[[#]]] &]];

TermToList[a_Times] := Join @@ (TermToList /@ List @@ a);
TermToList[Power[a_, n_Integer]] := ConstantArray[a, n] /; n > 0;
TermToList[a_] := {a};
SyntaxInformation@TermToList = {"ArgumentsPattern" -> {_}};

SeparateFactor[expr_, fns_List] := With[{
    l = FoldPairList[SeparateFactor, expr, fns, Identity]
}, Append[l[[All, 1]], l[[-1, 2]]]];
SeparateFactor[expr_, fn_] := Times @@@ Lookup[GroupBy[TermToList@expr, fn], {True, False}, 1];
SyntaxInformation@SeparateFactor = {"ArgumentsPattern" -> {_, _}};

GroupPolynomialBy[expr_List, fn_] := GroupPolynomialBy[#, fn] & /@ expr;
GroupPolynomialBy[expr_Plus, fn_] := KeySort@Merge[GroupPolynomialBy[#, fn] & /@ List @@ expr, Total];
GroupPolynomialBy[expr_, fn_] := Association[Rule @@ SeparateFactor[expr, fn]];
GroupPolynomialBy[fn_][expr_] := GroupPolynomialBy[expr, fn];
SyntaxInformation@GroupPolynomialBy = {"ArgumentsPattern" -> {_, _.}};

TermsToVecTable[terms_] := Association@Thread[terms -> IdentityMatrix@Length@terms];

PolynomialToVec[0, _, terms_] := ConstantArray[0, Length@terms];
PolynomialToVec[expr_, fn_, terms_] := GroupPolynomialBy[expr, fn] // KeyMap@TermsToVecTable@terms // KeyValueMap@Times // Total;
SyntaxInformation@PolynomialToVec = {"ArgumentsPattern" -> {_, _, _}};

CollectBy[expr_List, fn_, action_] := CollectBy[#, fn, action] & /@ expr;
CollectBy[expr_, fn_, action_] := GroupPolynomialBy[Expand@expr, fn] // Map@action // KeyValueMap@Times // Total;
CollectBy[fn_, action_][expr_] := CollectBy[expr, fn, action];
CollectBy[fn_][expr_] := CollectBy[expr, fn, Identity];
SyntaxInformation@CollectBy = {"ArgumentsPattern" -> {_, _., _.}};

ToStringInputForm[expr_String] := expr;
ToStringInputForm[expr_] := ToString@Unevaluated@InputForm@expr;
SetAttributes[ToStringInputForm, HoldAll];

AssignOpOfLHS[_SyntaxInformation] = "=";
AssignOpOfLHS[_] = ":=";
SetAttributes[AssignOpOfLHS, HoldAll];
WriteOneDownRule[stream_, Verbatim[HoldPattern][lhs_] :> rhs_] := WriteLine[
    stream,
    ToStringInputForm@lhs <> " " <> AssignOpOfLHS@lhs <> " " <> ToStringInputForm@rhs <> ";"
];
WriteOneUpRule[stream_, tag_, Verbatim[HoldPattern][lhs_] :> rhs_] := WriteLine[
    stream,
    ToStringInputForm@tag <> " /: " <>  ToStringInputForm@lhs <> " := " <> ToStringInputForm@rhs <> ";"
];
SetAttributes[{WriteOneDownRule, WriteOneUpRule}, HoldAll];

WriteDefinitionsAsDelayedValueTo[stream_, symbol_] := (
    Scan[WriteOneDownRule[stream, #] &, DownValues@symbol];
    Scan[WriteOneDownRule[stream, #] &, SubValues@symbol];
    Scan[WriteOneUpRule[stream, symbol, #] &, UpValues@symbol];
    Scan[WriteOneDownRule[stream, #] &, OwnValues@symbol];
    With[{attr = DeleteCases[Attributes@symbol, Temporary]}, If[Length@attr > 0,
        WriteLine[stream, "Attributes[" <> ToStringInputForm@symbol <> "] = " <> ToString@InputForm@attr <> ";"];
    ]];
    Scan[WriteOneDownRule[stream, #] &, DefaultValues@symbol];
    Scan[WriteOneUpRule[stream, symbol, #] &, FormatValues@symbol];
);
WriteDefinitionsAsDelayedValueTo[filename_String, expr_] := With[{outFile = OpenWrite@filename}, WithCleanup[
    WriteDefinitionsAsDelayedValueTo[outFile, expr],
    Close@outFile
]];
WriteDefinitionsAsDelayedValueTo[stream_, l_List] := WriteDefinitionsAsDelayedValueTo[stream, #] & /@ l;
SyntaxInformation@WriteDefinitionsAsDelayedValueTo = {"ArgumentsPattern" -> {_, _}};

SavableObjQ[expr_String] := ToExpression[expr, InputForm, SavableObjQ];
SavableObjQ[_] = True;
SetAttributes[SavableObjQ, HoldAll];
SyntaxInformation@SavableObjQ = {"ArgumentsPattern" -> {_}};

SaveNotebookData[] := SaveNotebookData@Select[Names["Global`*"], SavableObjQ];
SaveNotebookData[names_] := (
    SetDirectory[NotebookDirectory[]];
    With[{
        fname = FileBaseName[NotebookFileName[]] <> ".data.wl"
    },
        Quiet[DeleteFile[fname <> ".tmp"], DeleteFile::fdnfnd];
        WriteDefinitionsAsDelayedValueTo[fname <> ".tmp", names];
        Quiet[DeleteFile[fname], DeleteFile::fdnfnd];
        RenameFile[fname <> ".tmp", fname];
    ];
    ResetDirectory[];
);
SyntaxInformation@SaveNotebookData = {"ArgumentsPattern" -> {_.}};

LoadNotebookData[] := (
    SetDirectory[NotebookDirectory[]];
    Get[FileBaseName[NotebookFileName[]] <> ".data.wl"];
    ResetDirectory[];
);
SyntaxInformation[LoadNotebookData] = {"ArgumentsPattern" -> {}};

End[];

Protect @@ Names[$Context <> "*"];

EndPackage[];