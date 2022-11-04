result = TestReport[TestReport /@ {
    "xTension.wlt",
    "xDecomp.wlt"
}];

If[result["AllTestsSucceeded"],
    Print["All tests passed."]
,
    failed = result["TestsFailed"];
    Print[failed];
];
