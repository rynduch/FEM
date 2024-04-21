import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class Grid extends GlobalData {
  Node[] nodes = new Node[nodesNumber];
  Element[] elements = new Element[elementsNumber];

  public Grid() {
  }

  public Grid(GlobalData gd, Node[] n, Element[] e) {
    super(gd);
    nodes = n;
    elements = e;
  }

  public void count(UnivElement UE) {
    countH(UE);
    countHbc(UE);
    countP(UE);
    countC(UE);
    countHCdT();
    countCdTT0P();
    countCdT();
  }

  public void displayData() {
    System.out.println();
    System.out.println("Simulation Time: " + simulationTime);
    System.out.println("Simulation Step Time: " + simulationStepTime);
    System.out.println("Conductivity: " + conductivity);
    System.out.println("Alfa: " + alfa);
    System.out.println("Tot: " + tot);
    System.out.println("Initial Temp: " + initialTemp);
    System.out.println("Density: " + density);
    System.out.println("Specific Heat: " + specificHeat);
    System.out.println("Nodes number: " + nodesNumber);
    System.out.println("Elements number: " + elementsNumber);
    System.out.println("*Node");
    for (int i = 0; i < nodesNumber; i++) {
      System.out.println((i + 1) + ", " + nodes[i].x + ", " + nodes[i].y);
    }
    System.out.println("*Element");
    for (int i = 0; i < elementsNumber; i++) {
      System.out.print((i + 1));
      for (int j = 0; j < 4; j++) {
        System.out.print(", " + elements[i].ID[j]);
      }
      System.out.println();
    }
    System.out.println("*BC");
    for (int i = 0; i < nodesNumber; i++) {
      if (nodes[i].BC == 1) {
        System.out.print((i + 1) + "  ");
      }
    }
    System.out.println();
  }

  public void importData(String plik) {
    File file = new File(plik);
    try {
      Scanner scanner = new Scanner(file);
      while (scanner.hasNextLine()) {
        String line = scanner.nextLine();
        String[] parts = line.split(" ");
        switch (parts[0]) {
          case "SimulationTime" -> this.simulationTime = Integer.parseInt(parts[1]);
          case "SimulationStepTime" -> this.simulationStepTime = Integer.parseInt(parts[1]);
          case "Conductivity" -> this.conductivity = Integer.parseInt(parts[1]);
          case "Alfa" -> this.alfa = Integer.parseInt(parts[1]);
          case "Tot" -> this.tot = Integer.parseInt(parts[1]);
          case "InitialTemp" -> this.initialTemp = Integer.parseInt(parts[1]);
          case "Density" -> this.density = Integer.parseInt(parts[1]);
          case "SpecificHeat" -> this.specificHeat = Integer.parseInt(parts[1]);
          case "Nodes" -> {
            this.nodesNumber = Integer.parseInt(parts[2]);
            this.nodes = new Node[nodesNumber];
          }
          case "Elements" -> {
            this.elementsNumber = Integer.parseInt(parts[2]);
            this.elements = new Element[elementsNumber];
          }
          case "*Node" -> {
            for (int i = 0; i < nodesNumber; i++) {
              line = scanner.nextLine();
              parts = line.split(",");
              this.nodes[i] = new Node(Double.parseDouble(parts[1]), Double.parseDouble(parts[2]));
            }
          }
          case "*Element," -> {
            for (int i = 0; i < elementsNumber; i++) {
              line = scanner.nextLine();
              parts = line.split(",");
              this.elements[i] = new Element(new int[]{Integer.parseInt(parts[1].trim()), Integer.parseInt(parts[2].trim()), Integer.parseInt(parts[3].trim()), Integer.parseInt(parts[4].trim())});
            }
          }
          case "*BC" -> {
            line = scanner.nextLine();
            parts = line.split(",");
            for (String part : parts) {
              if (!part.isEmpty()) {
                int n = Integer.parseInt(part.trim());
                nodes[n - 1].setBC(1);
              }
            }
          }
        }
      }
      scanner.close();
    } catch (FileNotFoundException e) {
      System.out.println("Error. File not found.");
      e.printStackTrace();
    }
  }

  public void countJacobian(UnivElement UE, int id, int p) {

    double[][] jacobian = new double[2][2];

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        jacobian[i][j] = 0;
      }
    }

    jacobian[0][0] = UE.dNdKsi[p - 1][0] * nodes[elements[id - 1].ID[0] - 1].x + UE.dNdKsi[p - 1][1] * nodes[elements[id - 1].ID[1] - 1].x + UE.dNdKsi[p - 1][2] * nodes[elements[id - 1].ID[2] - 1].x + UE.dNdKsi[p - 1][3] * nodes[elements[id - 1].ID[3] - 1].x;
    jacobian[0][1] = UE.dNdKsi[p - 1][0] * nodes[elements[id - 1].ID[0] - 1].y + UE.dNdKsi[p - 1][1] * nodes[elements[id - 1].ID[1] - 1].y + UE.dNdKsi[p - 1][2] * nodes[elements[id - 1].ID[2] - 1].y + UE.dNdKsi[p - 1][3] * nodes[elements[id - 1].ID[3] - 1].y;
    jacobian[1][0] = UE.dNdEta[p - 1][0] * nodes[elements[id - 1].ID[0] - 1].x + UE.dNdEta[p - 1][1] * nodes[elements[id - 1].ID[1] - 1].x + UE.dNdEta[p - 1][2] * nodes[elements[id - 1].ID[2] - 1].x + UE.dNdEta[p - 1][3] * nodes[elements[id - 1].ID[3] - 1].x;
    jacobian[1][1] = UE.dNdEta[p - 1][0] * nodes[elements[id - 1].ID[0] - 1].y + UE.dNdEta[p - 1][1] * nodes[elements[id - 1].ID[1] - 1].y + UE.dNdEta[p - 1][2] * nodes[elements[id - 1].ID[2] - 1].y + UE.dNdEta[p - 1][3] * nodes[elements[id - 1].ID[3] - 1].y;

    UE.detJacobian[p - 1] = (jacobian[0][0] * jacobian[1][1]) - (jacobian[0][1] * jacobian[1][0]);

    double[][] inverseJacobian = new double[2][2];
    inverseJacobian[0][0] = jacobian[1][1];
    inverseJacobian[0][1] = -jacobian[0][1];
    inverseJacobian[1][0] = -jacobian[1][0];
    inverseJacobian[1][1] = jacobian[0][0];

    inverseJacobian[0][0] *= (1.0 / UE.detJacobian[p - 1]);
    inverseJacobian[0][1] *= (1.0 / UE.detJacobian[p - 1]);
    inverseJacobian[1][0] *= (1.0 / UE.detJacobian[p - 1]);
    inverseJacobian[1][1] *= (1.0 / UE.detJacobian[p - 1]);

    for (int j = 0; j < 4; j++) {
      UE.dNdX[p - 1][j] = inverseJacobian[0][0] * UE.dNdKsi[p - 1][j] + inverseJacobian[0][1] * UE.dNdEta[p - 1][j];
      UE.dNdY[p - 1][j] = inverseJacobian[1][0] * UE.dNdKsi[p - 1][j] + inverseJacobian[1][1] * UE.dNdEta[p - 1][j];
    }

  }

  public void displayJacobian(UnivElement UE, double[][] jacobian, int p) {
    System.out.println("Punkt " + p + ". ");
    System.out.println("- Jakobian:");
    for (int i = 0; i < 2; i++) {
      System.out.println("\t" + jacobian[i][0] + "\t" + jacobian[i][1]);
    }
    System.out.println("- detJakobian: " + UE.detJacobian[p - 1]);
  }

  public void countH(UnivElement UE) {
    for (int e = 1; e <= elementsNumber; e++) {
      for (int p = 1; p <= UE.numberOfIntegrationPoints; p++) {
        countJacobian(UE, e, p);
      }
      elements[e - 1].setH(UE, conductivity);
    }
  }

  public void countHbc(UnivElement UE) {
    for (int e = 1; e <= elementsNumber; e++) {
      elements[e - 1].setHbc(UE, nodes, alfa);
    }
  }

  public void countP(UnivElement UE) {
    for (int e = 1; e <= elementsNumber; e++) {
      elements[e - 1].setP(UE, nodes, alfa, tot);
    }
  }

  public void countC(UnivElement UE) {
    for (int e = 1; e <= elementsNumber; e++) {
      for (int p = 1; p <= UE.numberOfIntegrationPoints; p++) {
        countJacobian(UE, e, p);
      }
      elements[e - 1].setC(UE, density, specificHeat);
    }
  }

  public void countHCdT() {
    for (int e = 1; e <= elementsNumber; e++) {
      elements[e - 1].setHCdT(simulationStepTime);
    }
  }

  public void countCdTT0P() {
    for (int e = 1; e <= elementsNumber; e++) {
      elements[e - 1].setCdTT0P(simulationStepTime, initialTemp);
    }
  }

  public void countCdT() {
    for (int e = 1; e <= elementsNumber; e++) {
      elements[e - 1].setCdT(simulationStepTime);
    }
  }

  public void display() {
    System.out.println();
    System.out.println("----------------------------------------------------------------------------------------------------");
    System.out.println("                                            ELEMENTS");
    System.out.println("----------------------------------------------------------------------------------------------------");
    for (int e = 1; e <= elementsNumber; e++) {
      System.out.println();
      System.out.println("### Element " + e + " ###");
      System.out.println();
      elements[e - 1].displayH();
      elements[e - 1].displayHbc();
      elements[e - 1].displayP();
      elements[e - 1].displayC();
      elements[e - 1].displayHCdT();
      elements[e - 1].displayCdTT0P();
    }
  }
}
