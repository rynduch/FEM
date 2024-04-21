public class Aggregation {
  double[][] H_G;
  double[][] HHbc_G;
  double[] P_G;
  double[][] HCdT_G;
  double[] CdTT0P_G;
  double[][] CdT_G;
  Grid aGrid;

  public Aggregation(Grid G) {
    setH_G(G);
    setP_G(G);
    setHHbc_G(G);
    setHCdT_G(G);
    setCdTT0P_G(G);
    setCdT_G(G);
    aGrid = G;
  }

  public void display() {
    System.out.println();
    System.out.println("----------------------------------------------------------------------------------------------------");
    System.out.println("                                           AGGREGATION");
    System.out.println("----------------------------------------------------------------------------------------------------");
    System.out.println();
    displayH_G();
    displayHHbc_G();
    displayP_G();
    displayHCdT_G();
    displayCdTT0P_G();
    displayCdT_G();
    System.out.println();
  }

  public void setH_G(Grid G) {
    H_G = new double[G.nodesNumber][G.nodesNumber];
    for (int w = 0; w < G.elementsNumber; w++) {
      for (int k = 0; k < G.elementsNumber; k++) {
        H_G[w][k] = 0;
      }
    }
    for (int e = 0; e < G.elementsNumber; e++) {
      for (int w = 0; w < 4; w++) {
        for (int k = 0; k < 4; k++) {
          H_G[G.elements[e].ID[w] - 1][G.elements[e].ID[k] - 1] += G.elements[e].H[w][k];
        }
      }
    }
  }

  public void displayH_G() {
    System.out.println("Matrix [HG]: ");
    for (double[] doubles : H_G) {
      for (int w = 0; w < H_G.length; w++) {
        System.out.print("\t" + doubles[w]);
      }
      System.out.println();
    }
  }

  public void setHHbc_G(Grid G) {
    HHbc_G = new double[G.nodesNumber][G.nodesNumber];
    for (int w = 0; w < G.nodesNumber; w++) {
      System.arraycopy(H_G[w], 0, HHbc_G[w], 0, G.nodesNumber);
    }
    for (int e = 0; e < G.elementsNumber; e++) {
      for (int w = 0; w < 4; w++) {
        for (int k = 0; k < 4; k++) {
          HHbc_G[G.elements[e].ID[w] - 1][G.elements[e].ID[k] - 1] += G.elements[e].Hbc[w][k];
        }
      }
    }
  }

  public void displayHHbc_G() {
    System.out.println("Matrix [HHbcG]: ");
    for (double[] doubles : HHbc_G) {
      for (int w = 0; w < HHbc_G.length; w++) {
        System.out.print("\t" + doubles[w]);
      }
      System.out.println();
    }
  }

  public void setP_G(Grid G) {
    P_G = new double[G.nodesNumber];
    for (int i = 0; i < G.nodesNumber; i++) {
      P_G[i] = 0;
    }
    for (int e = 0; e < G.elementsNumber; e++) {
      for (int i = 0; i < 4; i++) {
        P_G[G.elements[e].ID[i] - 1] += G.elements[e].P[i];
      }
    }
  }

  public void displayP_G() {
    System.out.println("Vector [PG]: ");
    for (int i = 0; i < H_G.length; i++) {
      System.out.print("\t" + this.P_G[i]);
    }
    System.out.println();
  }

  public void setHCdT_G(Grid G) {
    HCdT_G = new double[G.nodesNumber][G.nodesNumber];
    for (int w = 0; w < G.elementsNumber; w++) {
      for (int k = 0; k < G.elementsNumber; k++) {
        HCdT_G[w][k] = 0;
      }
    }
    for (int e = 0; e < G.elementsNumber; e++) {
      for (int w = 0; w < 4; w++) {
        for (int k = 0; k < 4; k++) {
          HCdT_G[G.elements[e].ID[w] - 1][G.elements[e].ID[k] - 1] += G.elements[e].HCdT[w][k];
        }
      }
    }
  }

  public void displayHCdT_G() {
    System.out.println("Matrix [HCdT_G]: ");
    for (double[] doubles : HCdT_G) {
      for (int w = 0; w < HCdT_G.length; w++) {
        System.out.print("\t" + doubles[w]);
      }
      System.out.println();
    }
  }

  public void setCdTT0P_G(Grid G) {
    CdTT0P_G = new double[G.nodesNumber];
    for (int i = 0; i < G.nodesNumber; i++) {
      CdTT0P_G[i] = 0;
    }
    for (int e = 0; e < G.elementsNumber; e++) {
      for (int i = 0; i < 4; i++) {
        CdTT0P_G[G.elements[e].ID[i] - 1] += G.elements[e].CdTT0P[i];
      }
    }
  }

  public void setCdT_G(Grid G) {
    CdT_G = new double[G.nodesNumber][G.nodesNumber];
    for (int w = 0; w < G.nodesNumber; w++) {
      for (int k = 0; k < G.nodesNumber; k++) {
        CdT_G[w][k] = 0;
      }
    }
    for (int e = 0; e < G.elementsNumber; e++) {
      for (int w = 0; w < 4; w++) {
        for (int k = 0; k < 4; k++) {
          CdT_G[G.elements[e].ID[w] - 1][G.elements[e].ID[k] - 1] += G.elements[e].CdT[w][k];
        }
      }
    }
  }

  public void displayCdT_G() {
    System.out.println("Matrix [CdT_G]: ");
    for (double[] doubles : CdT_G) {
      for (int w = 0; w < CdT_G.length; w++) {
        System.out.print("\t" + doubles[w]);
      }
      System.out.println();
    }
  }

  public void displayCdTT0P_G() {
    System.out.println("Vector [CdTT0P_G]: ");
    for (double v : CdTT0P_G) {
      System.out.print("\t" + v);
    }
    System.out.println();
  }
}
