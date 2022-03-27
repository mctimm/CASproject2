const fs = require("fs");
const csv = require("csvtojson");

function getBest() {
  var store = new Array();
  let newstore = [];
  [10, "1e-07", 1000].forEach((name) => {
    csv()
      .fromFile("12_germinal_second_weight=" + name + "/BestAgainst1.csv")
      .then((jsonObj) => {
        jsonObj.forEach((newobj) => {
          let stuff = store.push({
            timestep: newobj.time,
            weight: name,
            average_score:
              Object.keys(newobj)
                .filter((x) => x.includes("G"))
                .reduce((x, y) => {
                  return x + parseFloat(newobj["" + y]);
                }, 0) / 12,
          });
        });
        fs.writeFileSync(
          "summary_data_all_BestAgainst1.csv.json",
          JSON.stringify(store)
        );
      });
  });
}

getBest();
