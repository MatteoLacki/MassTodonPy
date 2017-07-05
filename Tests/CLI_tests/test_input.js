var fs = require('fs')

var config = {
        "fasta":                "RPKPQQFFGLM",
        "precursor_charge":     3,
        "modifications":        {'C11':{'H':1,'O':-1,'N':1}},
        "cut_off":              100,
        "mz_prec":               .05
};

console.log(config)

fs.writeFile(
        '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/CLI_tests/config.json',
	JSON.stringify(config),
	'utf8',
	function(err) {
    	if (err) throw err;
    	console.log('complete');
    }
);
