BeginPackage["xTools`Misc`"];

(Unprotect[#]; ClearAll[#];) & /@ Names@{$Context <> "*", $Context <> "Private`*"};

FVariation::usage = "FVariation[expr, order] computes the functional variation.";
ConstantFunctions::usage = "ConstantFunctions is an option of FVariation that specifies constant functions.";
ConstantxTensors::usage = "ConstantxTensors is an option of FVariation that specifies constant tensors from xAct.";

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

End[];

Protect @@ Names[$Context <> "*"];

EndPackage[];