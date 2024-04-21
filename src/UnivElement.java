public class UnivElement {
  double[][] dNdKsi;
  double[][] dNdEta;
  double[][] dNdX;
  double[][] dNdY;
  double[][] N;
  int numberOfIntegrationPoints;
  double[] detJacobian;

  // Gauss table values
  final double val_n1 = 1.0 / Math.sqrt(3.0);
  final double val_n2 = Math.sqrt(3.0 / 5.0);
  final double val_1_n3 = Math.sqrt(3.0 / 7.0 + 2.0 / 7.0 * Math.sqrt(6.0 / 5.0));
  final double val_2_n3 = Math.sqrt(3.0 / 7.0 - 2.0 / 7.0 * Math.sqrt(6.0 / 5.0));

  // Integration points
  double[] gaussKsi_n1 = {-val_n1, val_n1, -val_n1, val_n1};
  double[] gaussEta_n1 = {-val_n1, -val_n1, val_n1, val_n1};
  double[] gaussKsi_n2 = {-val_n2, 0, val_n2, -val_n2, 0, val_n2, -val_n2, 0, val_n2};
  double[] gaussEta_n2 = {-val_n2, -val_n2, -val_n2, 0, 0, 0, val_n2, val_n2, val_n2};
  double[] gaussKsi_n3 = {-val_1_n3, -val_2_n3, val_2_n3, val_1_n3, -val_1_n3, -val_2_n3, val_2_n3, val_1_n3, -val_1_n3, -val_2_n3, val_2_n3, val_1_n3, -val_1_n3, -val_2_n3, val_2_n3, val_1_n3};
  double[] gaussEta_n3 = {-val_1_n3, -val_1_n3, -val_1_n3, -val_1_n3, -val_2_n3, -val_2_n3, -val_2_n3, -val_2_n3, val_2_n3, val_2_n3, val_2_n3, val_2_n3, val_1_n3, val_1_n3, val_1_n3, val_1_n3};

  // Weights
  double[] wN1 = {1, 1, 1, 1};
  double w1_N2 = 5.0 / 9.0, w2_N2 = 8.0 / 9.0, w3_N2 = 5.0 / 9.0;
  double[] wN2 = {w1_N2 * w1_N2, w2_N2 * w1_N2, w3_N2 * w1_N2, w1_N2 * w2_N2, w2_N2 * w2_N2, w3_N2 * w2_N2, w1_N2 * w3_N2, w2_N2 * w3_N2, w3_N2 * w3_N2};
  double w1_N3 = (18.0 - Math.sqrt(30.0)) / 36, w2_N3 = (18.0 + Math.sqrt(30.0)) / 36, w3_N3 = (18.0 + Math.sqrt(30.0)) / 36, w4_N3 = (18.0 - Math.sqrt(30.0)) / 36;
  double[] wN3 = {w1_N3 * w1_N3, w2_N3 * w1_N3, w3_N3 * w1_N3, w4_N3 * w1_N3, w1_N3 * w2_N3, w2_N3 * w2_N3, w3_N3 * w2_N3, w4_N3 * w2_N3, w1_N3 * w3_N3, w2_N3 * w3_N3, w3_N3 * w3_N3, w4_N3 * w3_N3, w1_N3 * w4_N3, w2_N3 * w4_N3, w3_N3 * w4_N3, w4_N3 * w4_N3};

  // Boundary integration points
  double[][] pcKsiBC_n1 = {{-val_n1, val_n1}, {1, 1}, {val_n1, -val_n1}, {-1, -1}};
  double[][] pcEtaBC_n1 = {{-1, -1}, {-val_n1, val_n1}, {1, 1}, {val_n1, -val_n1}};
  double[][] pcKsiBC_n2 = {{-val_n2, 0, val_n2}, {1, 1, 1}, {val_n2, 0, -val_n2}, {-1, -1, -1}};
  double[][] pcEtaBC_n2 = {{-1, -1, -1}, {-val_n2, 0, val_n2}, {1, 1, 1}, {val_n2, 0, -val_n2}};
  double[][] pcKsiBC_n3 = {{-val_1_n3, -val_2_n3, val_2_n3, val_1_n3}, {1, 1, 1, 1}, {val_1_n3, val_2_n3, -val_2_n3, -val_1_n3}, {-1, -1, -1, -1}};
  double[][] pcEtaBC_n3 = {{-1, -1, -1, -1}, {-val_1_n3, -val_2_n3, val_2_n3, val_1_n3}, {1, 1, 1, 1}, {val_1_n3, val_2_n3, -val_2_n3, -val_1_n3}};

  public class Surface {
    double[][] N_surface;

    public Surface() {
    }

    public void setN_surface(int S) {
      int n = (int) Math.sqrt(numberOfIntegrationPoints);
      N_surface = new double[n][4];
      if (numberOfIntegrationPoints == 4) {
        for (int i = 0; i < n; i++) {
          N_surface[i][0] = 0.25 * (1 - pcKsiBC_n1[S][i]) * (1 - pcEtaBC_n1[S][i]);
          N_surface[i][1] = 0.25 * (1 + pcKsiBC_n1[S][i]) * (1 - pcEtaBC_n1[S][i]);
          N_surface[i][2] = 0.25 * (1 + pcKsiBC_n1[S][i]) * (1 + pcEtaBC_n1[S][i]);
          N_surface[i][3] = 0.25 * (1 - pcKsiBC_n1[S][i]) * (1 + pcEtaBC_n1[S][i]);
        }
      } else if (numberOfIntegrationPoints == 9) {
        for (int i = 0; i < n; i++) {
          N_surface[i][0] = 0.25 * (1 - pcKsiBC_n2[S][i]) * (1 - pcEtaBC_n2[S][i]);
          N_surface[i][1] = 0.25 * (1 + pcKsiBC_n2[S][i]) * (1 - pcEtaBC_n2[S][i]);
          N_surface[i][2] = 0.25 * (1 + pcKsiBC_n2[S][i]) * (1 + pcEtaBC_n2[S][i]);
          N_surface[i][3] = 0.25 * (1 - pcKsiBC_n2[S][i]) * (1 + pcEtaBC_n2[S][i]);
        }
      } else if (numberOfIntegrationPoints == 16) {
        for (int i = 0; i < n; i++) {
          N_surface[i][0] = 0.25 * (1 - pcKsiBC_n3[S][i]) * (1 - pcEtaBC_n3[S][i]);
          N_surface[i][1] = 0.25 * (1 + pcKsiBC_n3[S][i]) * (1 - pcEtaBC_n3[S][i]);
          N_surface[i][2] = 0.25 * (1 + pcKsiBC_n3[S][i]) * (1 + pcEtaBC_n3[S][i]);
          N_surface[i][3] = 0.25 * (1 - pcKsiBC_n3[S][i]) * (1 + pcEtaBC_n3[S][i]);
        }
      }
    }
  }

  public Surface[] surfaces = new Surface[4];

  public void setSurfaces() {
    for (int i = 0; i < 4; i++) {
      surfaces[i] = new Surface();
      surfaces[i].setN_surface(i);
    }
  }

  public UnivElement(int n) {

    if (n != 2 && n != 3 && n != 4)
      System.out.println("Incorrect number of integration points.");
    else {
      numberOfIntegrationPoints = n * n;
      dNdKsi = new double[numberOfIntegrationPoints][4];
      dNdEta = new double[numberOfIntegrationPoints][4];
      dNdX = new double[numberOfIntegrationPoints][4];
      dNdY = new double[numberOfIntegrationPoints][4];
      detJacobian = new double[numberOfIntegrationPoints];

      if (numberOfIntegrationPoints == 4) {
        for (int i = 0; i < numberOfIntegrationPoints; i++) {
          dNdKsi[i][0] = -0.25 * (1 - gaussEta_n1[i]);
          dNdKsi[i][1] = 0.25 * (1 - gaussEta_n1[i]);
          dNdKsi[i][2] = 0.25 * (1 + gaussEta_n1[i]);
          dNdKsi[i][3] = -0.25 * (1 + gaussEta_n1[i]);
          dNdEta[i][0] = -0.25 * (1 - gaussKsi_n1[i]);
          dNdEta[i][1] = -0.25 * (1 + gaussKsi_n1[i]);
          dNdEta[i][2] = 0.25 * (1 + gaussKsi_n1[i]);
          dNdEta[i][3] = 0.25 * (1 - gaussKsi_n1[i]);
        }
        setN();
      } else if (numberOfIntegrationPoints == 9) {
        for (int i = 0; i < numberOfIntegrationPoints; i++) {
          dNdKsi[i][0] = -0.25 * (1 - gaussEta_n2[i]);
          dNdKsi[i][1] = 0.25 * (1 - gaussEta_n2[i]);
          dNdKsi[i][2] = 0.25 * (1 + gaussEta_n2[i]);
          dNdKsi[i][3] = -0.25 * (1 + gaussEta_n2[i]);
          dNdEta[i][0] = -0.25 * (1 - gaussKsi_n2[i]);
          dNdEta[i][1] = -0.25 * (1 + gaussKsi_n2[i]);
          dNdEta[i][2] = 0.25 * (1 + gaussKsi_n2[i]);
          dNdEta[i][3] = 0.25 * (1 - gaussKsi_n2[i]);
        }
        setN();
      } else if (numberOfIntegrationPoints == 16) {
        for (int i = 0; i < numberOfIntegrationPoints; i++) {
          dNdKsi[i][0] = -0.25 * (1 - gaussEta_n3[i]);
          dNdKsi[i][1] = 0.25 * (1 - gaussEta_n3[i]);
          dNdKsi[i][2] = 0.25 * (1 + gaussEta_n3[i]);
          dNdKsi[i][3] = -0.25 * (1 + gaussEta_n3[i]);
          dNdEta[i][0] = -0.25 * (1 - gaussKsi_n3[i]);
          dNdEta[i][1] = -0.25 * (1 + gaussKsi_n3[i]);
          dNdEta[i][2] = 0.25 * (1 + gaussKsi_n3[i]);
          dNdEta[i][3] = 0.25 * (1 - gaussKsi_n3[i]);
        }
        setN();
      }
    }
  }

  public void display() {
    System.out.println("Universal Elements: ");
    System.out.println("- dNdKsi");
    for (int i = 0; i < numberOfIntegrationPoints; i++) {
      System.out.println((i + 1) + "." + "\t" + dNdKsi[i][0] + "\t" + dNdKsi[i][1] + "\t" + dNdKsi[i][2] + "\t" + dNdKsi[i][3]);
    }
    System.out.println("- dNdEta");
    for (int i = 0; i < numberOfIntegrationPoints; i++) {
      System.out.println((i + 1) + "\t" + dNdEta[i][0] + "\t" + dNdEta[i][1] + "\t" + dNdEta[i][2] + "\t" + dNdEta[i][3]);
    }
    System.out.println("- N");
    for (int i = 0; i < numberOfIntegrationPoints; i++) {
      System.out.println((i + 1) + "\t" + N[i][0] + "\t" + N[i][1] + "\t" + N[i][2] + "\t" + N[i][3]);
    }
  }

  public void setN() {
    N = new double[numberOfIntegrationPoints][4];
    if (numberOfIntegrationPoints == 4) {
      for (int i = 0; i < numberOfIntegrationPoints; i++) {
        N[i][0] = 0.25 * (1 - gaussKsi_n1[i]) * (1 - gaussEta_n1[i]);
        N[i][1] = 0.25 * (1 + gaussKsi_n1[i]) * (1 - gaussEta_n1[i]);
        N[i][2] = 0.25 * (1 + gaussKsi_n1[i]) * (1 + gaussEta_n1[i]);
        N[i][3] = 0.25 * (1 - gaussKsi_n1[i]) * (1 + gaussEta_n1[i]);
      }
    } else if (numberOfIntegrationPoints == 9) {
      for (int i = 0; i < numberOfIntegrationPoints; i++) {
        N[i][0] = 0.25 * (1 - gaussKsi_n2[i]) * (1 - gaussEta_n2[i]);
        N[i][1] = 0.25 * (1 + gaussKsi_n2[i]) * (1 - gaussEta_n2[i]);
        N[i][2] = 0.25 * (1 + gaussKsi_n2[i]) * (1 + gaussEta_n2[i]);
        N[i][3] = 0.25 * (1 - gaussKsi_n2[i]) * (1 + gaussEta_n2[i]);
      }
    } else if (numberOfIntegrationPoints == 16) {
      for (int i = 0; i < numberOfIntegrationPoints; i++) {
        N[i][0] = 0.25 * (1 - gaussKsi_n3[i]) * (1 - gaussEta_n3[i]);
        N[i][1] = 0.25 * (1 + gaussKsi_n3[i]) * (1 - gaussEta_n3[i]);
        N[i][2] = 0.25 * (1 + gaussKsi_n3[i]) * (1 + gaussEta_n3[i]);
        N[i][3] = 0.25 * (1 - gaussKsi_n3[i]) * (1 + gaussEta_n3[i]);
      }
    }
  }
}


