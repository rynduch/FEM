public class GlobalData {
  int simulationTime;
  int simulationStepTime;
  int conductivity;
  int alfa;
  int tot;
  int initialTemp;
  int density;
  int specificHeat;
  int nodesNumber;
  int elementsNumber;

  public GlobalData() {
    this.simulationTime = 0;
    this.simulationStepTime = 0;
    this.conductivity = 0;
    this.alfa = 0;
    this.tot = 0;
    this.initialTemp = 0;
    this.density = 0;
    this.specificHeat = 0;
    this.nodesNumber = 4;
    this.elementsNumber = 1;
  }

  public GlobalData(GlobalData g) {
    simulationTime = g.simulationTime;
    simulationStepTime = g.simulationStepTime;
    conductivity = g.conductivity;
    alfa = g.alfa;
    tot = g.tot;
    initialTemp = g.initialTemp;
    density = g.density;
    specificHeat = g.specificHeat;
    nodesNumber = g.nodesNumber;
    elementsNumber = g.elementsNumber;
  }

  public GlobalData(int st, int sst, int c, int a, int t, int it, int d, int sh, int nn, int en) {
    this.simulationTime = st;
    this.simulationStepTime = sst;
    this.conductivity = c;
    this.alfa = a;
    this.tot = t;
    this.initialTemp = it;
    this.density = d;
    this.specificHeat = sh;
    this.nodesNumber = nn;
    this.elementsNumber = en;
  }
}
