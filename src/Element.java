import static java.lang.Math.sqrt;
import static java.lang.Math.pow;

public class Element {
  int[] ID;
  double[][] H = new double[4][4];
  double[][] Hbc = new double[4][4];
  double[] P = new double[4];
  double[][] C = new double[4][4];
  double[][] CdT = new double[4][4];
  double[][] HCdT = new double[4][4];
  double[] CdTT0P = new double[4];

  public Element(int[] id) {
    this.ID = id;
  }

  public static double[][] matrixHpc(UnivElement UE, int conductivity, int p) {

    double[][] Mx = new double[4][4];
    double[][] My = new double[4][4];

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        Mx[i][j] = UE.dNdX[p - 1][i] * UE.dNdX[p - 1][j];
        My[i][j] = UE.dNdY[p - 1][i] * UE.dNdY[p - 1][j];
      }
    }

    double[][] Hpc = new double[4][4];
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        Hpc[i][j] = (Mx[i][j] + My[i][j]) * (conductivity * UE.detJacobian[p - 1]);
      }
    }

    return Hpc;
  }

  public void setH(UnivElement UE, int conductivity) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        H[i][j] = 0;
      }
    }
    for (int p = 1; p <= UE.numberOfIntegrationPoints; p++) {
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          if (UE.numberOfIntegrationPoints == 4) {
            H[i][j] += (matrixHpc(UE, conductivity, p)[i][j] * UE.wN1[p - 1]);
          } else if (UE.numberOfIntegrationPoints == 9) {
            H[i][j] += (matrixHpc(UE, conductivity, p)[i][j] * UE.wN2[p - 1]);
          } else if (UE.numberOfIntegrationPoints == 16) {
            H[i][j] += (matrixHpc(UE, conductivity, p)[i][j] * UE.wN3[p - 1]);
          }
        }
      }
    }
  }

  public void displayH() {
    System.out.println("Matrix [H]: ");
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        System.out.print("\t" + H[i][j]);
      }
      System.out.println();
    }
  }

  public void setHbc(UnivElement UE, Node[] nodes, int alfa) {
    for (int k = 0; k < 4; k++) {
      for (int w = 0; w < 4; w++) {
        Hbc[k][w] = 0;
      }
    }
    double[] detJacobian = new double[]{0, 0, 0, 0};
    int numberOfBC = 0;
    int[] elementSurface = new int[]{0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
      if (nodes[this.ID[i] - 1].BC == 1 && nodes[this.ID[(i + 1) % 4] - 1].BC == 1) {
        numberOfBC++;
        elementSurface[i] = 1;
        detJacobian[i] = sqrt(pow(nodes[this.ID[(i + 1) % 4] - 1].x - nodes[this.ID[i] - 1].x, 2) + pow(nodes[this.ID[(i + 1) % 4] - 1].y - nodes[this.ID[i] - 1].y, 2)) / 2;
      }
    }
    if (numberOfBC != 0) {
      UE.setSurfaces();
      for (int S = 0; S < 4; S++) {
        if (elementSurface[S] == 1) {
          if (UE.numberOfIntegrationPoints == 4) {
            for (int k = 0; k < 4; k++) {
              for (int w = 0; w < 4; w++) {
                Hbc[k][w] += UE.surfaces[S].N_surface[0][w] * UE.surfaces[S].N_surface[0][k] * alfa * detJacobian[S];
                Hbc[k][w] += UE.surfaces[S].N_surface[1][w] * UE.surfaces[S].N_surface[1][k] * alfa * detJacobian[S];
              }
            }
          } else if (UE.numberOfIntegrationPoints == 9) {
            for (int k = 0; k < 4; k++) {
              for (int w = 0; w < 4; w++) {
                Hbc[k][w] += UE.surfaces[S].N_surface[0][w] * UE.surfaces[S].N_surface[0][k] * alfa * UE.w1_N2 * detJacobian[S];
                Hbc[k][w] += UE.surfaces[S].N_surface[1][w] * UE.surfaces[S].N_surface[1][k] * alfa * UE.w2_N2 * detJacobian[S];
                Hbc[k][w] += UE.surfaces[S].N_surface[2][w] * UE.surfaces[S].N_surface[2][k] * alfa * UE.w3_N2 * detJacobian[S];
              }
            }
          } else if (UE.numberOfIntegrationPoints == 16) {
            for (int k = 0; k < 4; k++) {
              for (int w = 0; w < 4; w++) {
                Hbc[k][w] += UE.surfaces[S].N_surface[0][w] * UE.surfaces[S].N_surface[0][k] * alfa * UE.w1_N3 * detJacobian[S];
                Hbc[k][w] += UE.surfaces[S].N_surface[1][w] * UE.surfaces[S].N_surface[1][k] * alfa * UE.w2_N3 * detJacobian[S];
                Hbc[k][w] += UE.surfaces[S].N_surface[2][w] * UE.surfaces[S].N_surface[2][k] * alfa * UE.w3_N3 * detJacobian[S];
                Hbc[k][w] += UE.surfaces[S].N_surface[3][w] * UE.surfaces[S].N_surface[3][k] * alfa * UE.w4_N3 * detJacobian[S];
              }
            }
          }
        }
      }
    }
  }

  public void displayHbc() {
    System.out.println("Matrix [Hbc]: ");
    for (int w = 0; w < 4; w++) {
      System.out.println("\t" + this.Hbc[w][0] + "\t" + this.Hbc[w][1] + "\t" + this.Hbc[w][2] + "\t" + this.Hbc[w][3]);
    }
  }

  public void setP(UnivElement UE, Node[] nodes, int alfa, int tot) {
    for (int i = 0; i < 4; i++) {
      P[i] = 0;
    }
    int numberOfBC = 0;
    int[] elementSurface = new int[]{0, 0, 0, 0};
    double[] detJacobian = new double[]{0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
      if (nodes[this.ID[i] - 1].BC == 1 && nodes[this.ID[(i + 1) % 4] - 1].BC == 1) {
        numberOfBC++;
        elementSurface[i] = 1;
        detJacobian[i] = sqrt(pow(nodes[this.ID[(i + 1) % 4] - 1].x - nodes[this.ID[i] - 1].x, 2) + pow(nodes[this.ID[(i + 1) % 4] - 1].y - nodes[this.ID[i] - 1].y, 2)) / 2;
      }
    }
    if (numberOfBC != 0) {
      for (int S = 0; S < 4; S++) {
        if (elementSurface[S] == 1) {
          UE.setSurfaces();
          if (UE.numberOfIntegrationPoints == 4) {
            for (int i = 0; i < 4; i++) {
              P[i] += UE.surfaces[S].N_surface[0][i] * alfa * tot * detJacobian[S];
              P[i] += UE.surfaces[S].N_surface[1][i] * alfa * tot * detJacobian[S];

            }
          } else if (UE.numberOfIntegrationPoints == 9) {
            for (int i = 0; i < 4; i++) {
              P[i] += UE.surfaces[S].N_surface[0][i] * alfa * tot * UE.w1_N2 * detJacobian[S];
              P[i] += UE.surfaces[S].N_surface[1][i] * alfa * tot * UE.w2_N2 * detJacobian[S];
              P[i] += UE.surfaces[S].N_surface[2][i] * alfa * tot * UE.w3_N2 * detJacobian[S];

            }
          } else if (UE.numberOfIntegrationPoints == 16) {
            for (int i = 0; i < 4; i++) {
              P[i] += UE.surfaces[S].N_surface[0][i] * alfa * tot * UE.w1_N3 * detJacobian[S];
              P[i] += UE.surfaces[S].N_surface[1][i] * alfa * tot * UE.w2_N3 * detJacobian[S];
              P[i] += UE.surfaces[S].N_surface[2][i] * alfa * tot * UE.w3_N3 * detJacobian[S];
              P[i] += UE.surfaces[S].N_surface[3][i] * alfa * tot * UE.w4_N3 * detJacobian[S];
            }
          }
        }
      }
    }
  }

  public void displayP() {
    System.out.println("Vector [P]: ");
    System.out.println("\t" + this.P[0] + "\t" + this.P[1] + "\t" + this.P[2] + "\t" + this.P[3]);
  }

  public static double[][] matrixCpc(UnivElement UE, int density, int specificHeat, int p) {

    double[][] M = new double[4][4];

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        M[i][j] = UE.N[p - 1][i] * UE.N[p - 1][j];
      }
    }

    double[][] Cpc = new double[4][4];

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        Cpc[i][j] = M[i][j] * density * specificHeat * UE.detJacobian[p - 1];
      }
    }
    return Cpc;
  }

  public void setC(UnivElement UE, int density, int specificHeat) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        C[i][j] = 0;
      }
    }
    for (int p = 1; p <= UE.numberOfIntegrationPoints; p++) {
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          if (UE.numberOfIntegrationPoints == 4) {
            C[i][j] += (matrixCpc(UE, density, specificHeat, p)[i][j] * UE.wN1[p - 1]);
          } else if (UE.numberOfIntegrationPoints == 9) {
            C[i][j] += (matrixCpc(UE, density, specificHeat, p)[i][j] * UE.wN2[p - 1]);
          } else if (UE.numberOfIntegrationPoints == 16) {
            C[i][j] += (matrixCpc(UE, density, specificHeat, p)[i][j] * UE.wN3[p - 1]);
          }
        }
      }
    }
  }

  public void displayC() {
    System.out.println("Matrix [C]: ");
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        System.out.print("\t" + C[i][j]);
      }
      System.out.println();
    }
  }

  public void displayCpc(double[][] Cpc, int p) {
    System.out.println("Matrix [Cpc] " + p + ":");
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        System.out.print("\t" + Cpc[i][j]);
      }
      System.out.println();
    }
  }

  public void setHCdT(int simulationStepTime) {
    for (int w = 0; w < 4; w++) {
      for (int k = 0; k < 4; k++) {
        HCdT[w][k] = (C[w][k] / simulationStepTime) + H[w][k] + Hbc[w][k];
      }
    }
  }

  public void displayHCdT() {
    System.out.println("Matrix [HCdT]:");
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        System.out.print("\t" + HCdT[i][j]);
      }
      System.out.println();
    }
  }

  public void setCdTT0P(int simulationStepTime, int initialTemp) {
    for (int w = 0; w < 4; w++) {
      CdTT0P[w] = 0;
    }
    for (int w = 0; w < 4; w++) {
      for (int k = 0; k < 4; k++) {
        CdTT0P[w] += (C[w][k] / simulationStepTime) * initialTemp;
      }
      CdTT0P[w] += P[w];
    }
  }

  public void displayCdTT0P() {
    System.out.println("Vector [CdTT0P]: ");
    System.out.print("t1");
    for (int j = 0; j < 4; j++) {
      System.out.print("\t" + CdTT0P[j]);
    }
    System.out.println();
  }

  public void setCdT(int simulationStepTime) {
    for (int w = 0; w < 4; w++) {
      for (int k = 0; k < 4; k++) {
        CdT[w][k] = C[w][k] / simulationStepTime;
      }
    }
  }
}
