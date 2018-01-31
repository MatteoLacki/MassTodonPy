json_of_a_bitch = JSON.stringify([10, 11, 'adav'])

var fs = require('fs');

fs.writeFile("test.json", json_of_a_bitch, function(err) {
    if(err) {
        return console.log(err);
    }
});
