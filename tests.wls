result = TestReport["xTension.wlt"];

If[result["AllTestsSucceeded"],
    Print["All tests passed."],
    failed = result["TestsFailed"];
    Print[failed];
];
