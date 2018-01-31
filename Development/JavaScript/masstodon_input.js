input = {
    'spectrum': '/Users/matteo/Documents/MassTodon/Tests/ReadingJson/spectrum.txt', // path to the spectrum file - required
    'output': '/Users/matteo/Documents/MassTodon/Tests/ReadingJson/output/',        // path to output - required
    'fasta': 'RPKPQQFFGLM', // FASTA sequence - required
    'charge': 3,            // charge - required
    'mz_tol': 0.05,         // tolerance in m/z ratios - required
    'name': 'substance P',  // name - optional
    'modifications': {      // modifications - optional
        10: {
            'C_carbo': {'H': 1, 'N': 1, 'O':-1} 
        }
    },
    'fragmentation_type': 'cz',         // type of fragmentation: optional (right now supports only cz!!!)
    'blocked_fragments': 'c0z10',       // blocked fragments: optional 
    'distance_charges': 4,              // distance between charges: optional
    'min_intensity': 100.0,             // cut-off intensity threshold: optional
    'percent_top_peaks': 0.999,         // relative cut-off for joint intensity: optional
    'deconvolution_method': 'Matteo',   // to choose from 'Matteo' and 'Ciacho_Wanda': optional - default to Matteo
    'joint_probability': 0.999,         // joint probability of the smallest isotopic distribution: optional
    'min_prob_per_molecule': 0.7,       // criterion for a molecule to be included in the deconvolution graph: optional
    '_max_buffer_len': 0.5,             // extremely optional : extremely optional
    '_L1_flow': 0.001,                  // L1 penalty on flows : extremely optional
    '_L2_flow': 0.001,                  // L2 penalty on flows : extremely optional
    '_L1_intensity': 0.001,             // L1 penalty on intensity : extremely optional
    '_L2_intensity': 0.001,             // L2 penalty on intensity : extremely optional
    '_max_times': 10,                   // how many times to rerun CVXOPT : extremely optional
    '_maxiters': 1000,                  // maximal number of iterations of the iterative CVXOPT optimization routing : extremely optional
    'sigma2': 0.1,                      // the variance of the experimental peak's m/z ratio : optional
    'ni2': 0.1                          // the variance of the theoretic isotopologue's m/z ratio: optional
}

console.log(input)

json_of_a_bitch = JSON.stringify(input)

var fs = require('fs');

fs.writeFile("masstodon_input.json", json_of_a_bitch, function(err) {
    if (err) {
        return console.log(err);
    }
});