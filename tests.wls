xTools`$SkipLongTest = !MemberQ[$ScriptCommandLine[[2 ;;]], "--all"];

result = TestReport[TestReport /@ FileNames[__ ~~ ".wlt", DirectoryName@$InputFileName]];

If[result["AllTestsSucceeded"],
    Print["All tests passed."]
,
    failed = result["TestsFailed"];
    Print[failed /. Removed[name_] :> name];
];
