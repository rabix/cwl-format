key_order_dict = {

    "generic-ordering": [
        "id",
        "label",
        "name",
        "doc",
        "class",
        "type",
        "format",
        "default",
        "secondaryFiles",
        "inputBinding",
        "prefix",
        "position",
        "valueFrom",
        "separate",
        "itemSeparator",
        "shellQuote",
        "outputBinding",
        "glob",
        "outputEval",
        "loadContents",
        "loadListing",
        "dockerPull",

        # WorkflowStep
        "in",
        "scatter",
        "scatterMethod",
        "run",
        "when",
        "out",
        "requirements",
        "hints",

        # WorkflowStepInput
        "source",
        "outputSource"
        "linkMerge"
    ],

    "CommandLineTool": [
        "class",
        "cwlVersion",
        "label",
        "doc",
        "requirements",
        "inputs",
        "outputs",
        "baseCommand",
        "arguments",
        "stdout",
        "stderr",
        "hints",
        "id"
    ],

    "ExpressionTool": [
        "class",
        "cwlVersion",
        "label",
        "doc",
        "requirements",
        "inputs",
        "outputs",
        "expression",
        "hints",
        "id"
    ],

    "Workflow": [
        "class",
        "cwlVersion",
        "label",
        "doc",
        "requirements",
        "inputs",
        "outputs",
        "steps",
        "hints",
        "id"
    ]
}
