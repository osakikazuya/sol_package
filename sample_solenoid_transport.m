(* ::Package:: *)

(*RungeKutta Method 4th*)

 Off[NIntegrate::inumr];
 Off[NDSolveValue::ecboo];
 Off[NDSolveValue::nbnum1];

 $HistoryLength = 3;
 ClearAll[SOL, Xinit, Yinit, Pxinit, Pyinit];
 Needs["DifferentialEquations`NDSolveProblems`"];
 
 R = 0.08;(*Radius of solenoid coil[m]*)
 
 l = 0.3;(*Lenght of solenoid coil[m]*)
 
 Lz = 20;(*Simulation area of longitudinal direction[m]*)
 
 clight = 2.99792458*10^8;(*Speed of light[m/s]*)
 
 m = 238*1.66053892*10^-27;(*Mass of ion[kg]*)
 
 q = 1.60217657*10^-19;(*Elementary chrge*)
 
 Inum = 34;
 
 Ene = 3*10^6;(*Kinetic energy*)
 
 gamma = q*Ene/(m*clight^2) + 1;(*Lorentz factor*)
 
 v = clight*Sqrt[1 - 1/gamma^2];(*Velosity of ion*)
 
 Xinit[1] = {0.002, 0.004};(*Initial x coodinate[m]*)
 
 Yinit[1] = {0.002, 0};(*Initial y coodinate[m]*)
 
 Pxinit[1] = {100, 0};(*Initial vx[m/s]*)
 
 Pyinit[1] = {100, 100};(*Initial vy[m/s]*)
 
 BzCent[r_, z_] := NIntegrate[
                              1/(2 \[Pi]) ((R^2 - R*r*Cos[\[Theta]])/(
                                                                      R^2 + r^2 - 
                                                                      2 R*r*Cos[\[Theta]]))*((z + l/2)/
                                                                                             Sqrt[(z + l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]] - (z - l/2)/
                                                                                             Sqrt[(z - l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]]), {\[Theta], 0, 2 \[Pi]}];
 
 Bz0 = BzCent[0, 0];
 
 Br[r_, z_] := NIntegrate[
                          1/(2 \[Pi])*((R*Cos[\[Theta]])/
                                       Sqrt[(z - l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]] - (
                                                                                              R*Cos[\[Theta]])/
                                       Sqrt[(z + l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]])/
                          Bz0, {\[Theta], 0, 2 \[Pi]}];(*Br not scaled*)
 
 Bz[r_, z_] := NIntegrate[
                          1/(2 \[Pi]) ((R^2 - R*r*Cos[\[Theta]])/(
                                                                  R^2 + r^2 - 
                                                                  2 R*r*Cos[\[Theta]]))*((z + l/2)/
                                                                                         Sqrt[(z + l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]] - (z - l/2)/
                                                                                         Sqrt[(z - l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]])/
                          Bz0, {\[Theta], 0, 2 \[Pi]}];
 
 Vari = {x, y, z};(*Variable of equation*)
 
 Time = {t, 0, Lz/v};(*Simulation time*)
 
 Table[Eqs[i] = {m*gamma*x''[t] == Inum*q*(y'[t]*Bz[Sqrt[x[t]^2 + y[t]^2], z[t]] - z'[t]*y[t]/Sqrt[x[t]^2 + y[t]^2]*Br[Sqrt[x[t]^2 + y[t]^2], z[t]]), 
                 m*gamma*y''[t] == Inum*q*(z'[t]*x[t]/Sqrt[x[t]^2 + y[t]^2]*Br[Sqrt[x[t]^2 + y[t]^2], z[t]] - x'[t]*Bz[Sqrt[x[t]^2 + y[t]^2], z[t]]), 
                 m*gamma*z''[t] == Inum*q*(x'[t]*y[t]/Sqrt[x[t]^2 + y[t]^2]*Br[Sqrt[x[t]^2 + y[t]^2], z[t]] - y'[t] x[t]/Sqrt[x[t]^2 + y[t]^2]*
                                   Br[Sqrt[x[t]^2 + y[t]^2], z[t]]), x[0] == Xinit[1][[i]], 
                 x'[0] == Pxinit[1][[i]], y[0] == Yinit[1][[i]], y'[0] == Pyinit[1][[i]], z[0] == -Lz/2, z'[0] == v}, {i, 1, 2}];(*Equation of motion*)
                 Table[SOL[i] = NDSolveValue[Eqs[i], Vari, Time, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 4}, MaxSteps -> 10^5], {i, 1, 2}];
                 (*Solve equation of motion numerically*)

Table[xy[i] = Table[{SOL[i][[1]][t], SOL[i][[2]][t]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
Table[xz[i] = Table[{SOL[i][[3]][t], SOL[i][[1]][t]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
Table[yz[i] = Table[{SOL[i][[3]][t], SOL[i][[2]][t]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];

Print[Text[Style["Plot of Particle Orbit in Real Space", 16, Bold, Black]]];
    
Print[XY[1] = ListPlot[Table[xy[i], {i, 1, 2}], PlotRange -> {{-0.03, 0.03}, {-0.03, 0.03}}, FrameLabel -> {"x_coordinate[m]", "y_coordinate[m]"}, 
                       AspectRatio -> 1, LabelStyle -> Directive[15], PlotStyle -> PointSize[0.009], Joined -> True, Frame -> True]];
            
Print[XZ[1] = ListPlot[Table[xz[i], {i, 1, 2}], PlotRange -> {{-10, 10}, {-0.2, 0.2}}, FrameLabel -> {"z_coordinate[m]", "x_coordinate[m]"}, 
                LabelStyle -> Directive[15], PlotStyle -> PointSize[0.009], Joined -> True, Frame -> True]];
            
Print[YZ[1] = ListPlot[Table[yz[i], {i, 1, 2}], PlotRange -> {{-10, 10}, {-0.2, 0.2}}, FrameLabel -> {"z_coordinate[m]", "y_coordinate[m]"}, 
                 LabelStyle -> Directive[15], PlotStyle -> PointSize[0.009], Joined -> True, Frame -> True]];

Athe[r_, z_] := NIntegrate[(R^2*r)/(2 \[Pi]) ((Sin[\[Theta]]*Sin[\[Theta]])/(R^2 + r^2 - 2 R*r*Cos[\[Theta]]))*((z + l/2)/
                            Sqrt[(z + l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]] - (z - l/2)/Sqrt[(z - l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]])/(Bz0), {\[Theta], 0, 2 \[Pi]}];

Table[Invari[i] = Table[{SOL[i][[3]][t], m*gamma*(SOL[i][[1]][t]*SOL[i][[2]]'[t] - SOL[i][[2]][t]*SOL[i][[1]]'[t]) + 
                    Inum*q*Athe[Sqrt[SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2], SOL[i][[3]][t]]*Sqrt[SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2]}, {t,0, Lz/v, Lz/v/2000}], {i, 1, 2}];

Table[MAM[i] = Table[{SOL[i][[3]][t], m*gamma*(SOL[i][[1]][t]*SOL[i][[2]]'[t] - SOL[i][[2]][t]*SOL[i][[1]]'[t])}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
                    (*Mechanical_Anguler_momentum vs z*)
                    
Table[PT[i] = Table[{SOL[i][[3]][t], Inum*q*Athe[Sqrt[SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2], SOL[i][[3]][t]]*Sqrt[SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2]}, {t,0, Lz/v, Lz/v/2000}],{i, 1, 2}];

Table[Gamma3d[i] = Table[{SOL[i][[3]][t], 1/Sqrt[1 - (SOL[i][[1]]'[t]^2 + SOL[i][[2]]'[t]^2 + SOL[i][[3]]'[t]^2)/clight^2]}, {t, 0, Lz/v, Lz/v/2000}], {i,1, 2}];

Table[Beta3d[i] = Table[{SOL[i][[3]][t], Sqrt[(SOL[i][[1]]'[t]^2 + SOL[i][[2]]'[t]^2 + SOL[i][[3]]'[t]^2)]/clight}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];

Table[Gammaz[i] = Table[{SOL[i][[3]][t], 1/Sqrt[1 - (SOL[i][[3]]'[t]^2)/clight^2]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];

Table[Betaz[i] = Table[{SOL[i][[3]][t], SOL[i][[3]]'[t]/clight}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
                    (*Beta and gamma vs z*)
Table[Rorb[i] = Table[{SOL[i][[3]][t], Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
                    (*r vs z*)
Table[Rdorb[i] = Table[{SOL[i][[3]][t], ((SOL[i][[1]][t]*SOL[i][[1]]'[t] + SOL[i][[2]][t]*SOL[i][[2]]'[t])/Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)]*1/SOL[i][[3]]'[t])}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];

Print[Text[Style["Plot of P\[Theta], Kinetic term, and Potential term", 16, Bold, Black]]];
                                         
Print[Ps[1] = ListPlot[Table[Invari[i], {i, 1, 2}], PlotRange -> {{-5, 5}, {-0.3*10^-25, 5*10^-25}}, FrameLabel -> {"z_coordinate[m]", "\!\(\*SubscriptBox[\(P\), \(\[Theta]\)]\)"}, 
                       LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];
Print[Kt[1] = ListPlot[Table[MAM[i], {i, 1, 2}], PlotRange -> {{-5, 5}, {-2*10^-22, 2*10^-22}}, FrameLabel -> {"z_coordinate[m]", "Kitetic term"}, 
                       LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];
Print[Pt[1] = ListPlot[Table[PT[i], {i, 1, 2}], PlotRange -> {{-5, 5}, {-2*10^-22, 2*10^-22}}, FrameLabel -> {"z_coordinate[m]", "Potential term"}, 
                       LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];

Print[Text[Style["Plot of Lorentz factors", 16, Bold, Black]]];                                         

Print[LG[1] = Show[ListPlot[Table[Gamma3d[i], {i, 1, 2}], PlotRange -> {{-5, 5}, {1.000013525, 1.000013537}}, FrameTicks -> {{{{1.000013525, "1.000013525"}, {1.00001353, "1.00001353"}, {1.000013535, "1.000013535"}}, None}, {{-5, -2.5, 0, 2.5, 5}, None}}, FrameLabel -> {"z_coordinate[m]", "Lorenz factor \[Gamma]"}, LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}],
            ListPlot[Table[Gammaz[i], {i, 1, 2}],PlotRange -> {{-5, 5}, {1.0000135, 1.0000136}}, FrameTicks -> {{{{1.0000135, "1.0000135"}, {1.00001355, 
                    "1.00001355"}, {1.0000136, "1.0000136"}}, None}, {{-5, -2.5, 0, 2.5, 5}, None}}, FrameLabel -> {"z_coordinate[m]", "Lorenz factor \[Gamma]"}, 
                     LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]]];

Print[LB[1] = Show[ListPlot[Table[Beta3d[i], {i, 1, 2}], PlotRange -> {{-5, 5}, {0.005201, 0.005203}}, FrameTicks -> {{{{0.005201, "0.005201"},{0.005202, "0.005202"}, {0.005203, "0.005203"}},  
                    None}, {{-5, -2.5, 0, 2.5, 5}, None}}, FrameLabel -> {"z_coordinate[m]", "Lorenz factor \[Beta]"}, LabelStyle -> Directive[15], Frame -> True, Axes -> None, 
                    Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}],
             ListPlot[Table[Betaz[i], {i, 1, 2}], PlotRange -> {{-5, 5}, {0.00519, 0.00521}}, FrameLabel -> {"z_coordinate[m]", "Lorenz factor \!\(\*SubscriptBox[\(\[Beta]\), \(b\)]\)"}, 
                      LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]]];                                         
                    
Print[Text[Style["Plot of r and dr/dz orbit", 16, Bold, Black]]];
                                         
Print[Ro[1] = ListPlot[Table[Rorb[i], {i, 1, 2}], PlotRange -> {{-5, 5}, {0, 2*10^-2}}, FrameLabel -> {"z_coordinate[m]", "r_coordinate[m]"}, 
                       LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];

Print[Rdo[1] = ListPlot[Table[Rdorb[i], {i, 1, 2}], PlotRange -> {{-5, 5}, {-4*10^-2, 4*10^-2}}, FrameLabel -> {"z_coordinate[m]", "dr/dz"}, 
                        LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];                                         
                                    
Table[RKick[i, t_] := SOL[i][[3]]'[t]*((Inum*q*Bz[Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)], SOL[i][[3]][t]])/(m*gamma)*(SOL[i][[1]][t]*SOL[i][[2]]'[t] - 
                        SOL[i][[2]][t]*SOL[i][[1]]'[t])/(Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)]*SOL[i][[3]]'[t]^2) + 
                        1/((SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)^(3/2)*SOL[i][[3]]'[t]^2)*(-(SOL[i][[1]][t]*SOL[i][[1]]'[t] + 
                        SOL[i][[2]][t]*SOL[i][[2]]'[t])^2 + (SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)*(SOL[i][[1]]'[t]^2 + 
                        SOL[i][[2]]'[t]^2)) - (Inum*q*Br[Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)], SOL[i][[3]][t]])/(m*gamma*(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)*
                        SOL[i][[3]]'[t]^3)*(SOL[i][[1]][t]*SOL[i][[1]]'[t] + SOL[i][[2]][t]*SOL[i][[2]]'[t])*(SOL[i][[2]][t]*SOL[i][[1]]'[t] - SOL[i][[1]][t]*
                        SOL[i][[2]]'[t]) - ((m*gamma*(SOL[i][[1]][0]*SOL[i][[2]]'[0] - SOL[i][[2]][0]*SOL[i][[1]]'[0]) + 
                        Inum*q*Athe[Sqrt[SOL[i][[1]][0]^2 + SOL[i][[2]][0]^2], SOL[i][[3]][0]]*Sqrt[SOL[i][[1]][0]^2 + SOL[i][[2]][0]^2])/(m*gamma*v))^2/(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)^(3/2)), {i, 1,2}];
                                                                            
Table[RDD[i] = NIntegrate[RKick[i, t], {t, 0, Lz/v}], {i, 1, 2}];
      
ClearAll[r0];
Bz00[z_] = (z + l/2)/Sqrt[(z + l/2)^2 + R^2] - (z - l/2)/Sqrt[(z - l/2)^2 + R^2];(*Bz field*)
                                         
Bz1[z_] := ((z + l/2)/Sqrt[(z + l/2)^2 + R^2] - (z - l/2)/Sqrt[(z - l/2)^2 + R^2])/Bz00[0];(*Scaled Bz field*)
                                         
Bz2[z_] := Bz1'[z];(*Scaled Bz field*)
                                         
Bz3[z_] := Bz2'[z];(*Scaled Bz field*)
                                         
Bz4[z_] := Bz3'[z];(*Scaled Bz field*)
                                         
Bz5[z_] := Bz4'[z];(*Scaled Bz field*)
                                         
Bz6[z_] := Bz5'[z];(*Scaled Bz field*)

Brho = m*gamma*v/(Inum*q);
                                    
r0[1] = Table[Sqrt[Xinit[1][[i]]^2 + Yinit[1][[i]]^2], {i, 1, 2}];      
      
Table[RD1[i] = -1/4*NIntegrate[(Bz1[z]/Brho)^2*r0[1][[i]], {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];
                                         
Table[RD2[i] = -1/8*NIntegrate[(Bz2[z]/Brho)^2*r0[1][[i]]^3, {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];
                                         
Table[RD3[i] = -5/256*NIntegrate[(Bz3[z]/Brho)^2*r0[1][[i]]^5, {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];
                                         
Table[RD4[i] = -7/4608*NIntegrate[(Bz4[z]/Brho)^2*r0[1][[i]]^7, {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];
                                         
Table[RD5[i] = -7/98304*NIntegrate[(Bz5[z]/Brho)^2*r0[1][[i]]^9, {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];

Print[Text[Style["Comparison of Radial kick between the theory and the simulation", 16, Bold, Black]]];                                         
                                         
Print[Rdd[1] = Show[
              ListPlot[Table[{Invari[i][[1, 2]], (RD1[i] + RD2[i] + RD3[i] + RD4[i] + RD5[i])}, {i, 1, 2}], PlotRange -> {{-0.5*10^-25, 3.5*10^-25}, {-0.05, 0.05}}, 
                       FrameTicks -> {{{{-0.04, "-0.04"}, {-0.02, "-0.02"}, {0, "0"}, {0.02, "0.02"}, {0.04, "0.04"}}, None}, {{{0, "0"}, {1*10^-25, 
                        "1*\!\(\*SuperscriptBox[\(10\), \(-25\)]\)"}, {2*10^-25, "2*\!\(\*SuperscriptBox[\(10\), \(-25\)]\)"}, {3*10^-25, 
                        "3*\!\(\*SuperscriptBox[\(10\), \(-25\)]\)"}}, None}}, PlotStyle -> {Blue}, FrameLabel -> {"\!\(\*SubscriptBox[\(P\), \(\[Theta]\)]\)", 
                        "\[CapitalDelta]r'"}, LabelStyle -> Directive[15], Frame -> True,Axes -> None],
                        ListPlot[Table[{Invari[j][[1, 2]], RDD[j]}, {j, 1, 2}], PlotRange -> {{-0.5*10^-25, 3.5*10^-25}, {-0.05, 0.05}}, 
                        PlotStyle -> {Red}, FrameLabel -> {"\!\(\*SubscriptBox[\(P\), \(\[Theta]\)]\)","\[CapitalDelta]r'"}, LabelStyle -> Directive[15], Frame -> True,
                                 Axes -> None]]];
      
      
      
      
      

