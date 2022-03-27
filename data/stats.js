const fs = require("fs");
const csv = require("csvtojson");

function getBest() {
  var store = new Array();
  let newstore = [];
  for (let index = 0; index < 25; index++) {
    csv()
      .fromFile(`./12_germinal_${index}_transport/BestAgainst1.csv`)
      .then((jsonObj) => {
        jsonObj.forEach((newobj) => {
          let stuff = store.push({
            timestep: newobj.time,
            transport: index,
            average_score:
              Object.keys(newobj)
                .filter((x) => x.includes("G"))
                .reduce((x, y) => {
                  return x + parseFloat(newobj["" + y]);
                }, 0) / 12,
          });
        });
        fs.writeFileSync(
          "summary_data___BestAgainst1.json",
          JSON.stringify(store)
        );
      });
  }
}

getBest();
