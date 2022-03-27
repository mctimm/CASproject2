const fs = require("fs");
const csv = require("csvtojson");

function getBest() {
  var store = new Array();
  let newstore = [];
  for (let index = 2; index < 49; index += 2) {
    csv()
      .fromFile(`./12_germinal_${index}_virus/BestAgainst1.csv`)
      .then((jsonObj) => {
        jsonObj.forEach((newobj) => {
          let stuff = store.push({
            timestep: newobj.time,
            virus: index,
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
