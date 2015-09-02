(* ::Package:: *)

(*RungeKutta Method 4th*)
 Print[Text[Style["Single Particle Simulation through a Nonlinear Solenoid using 4th RungeKutta Method", 24, Bold, Black]]];
 Off[NIntegrate::inumr];
 Off[NDSolveValue::ecboo];
 Off[NDSolveValue::nbnum1];

 $HistoryLength = 3;
 ClearAll[SOL, Xinit, Yinit, Pxinit, Pyinit];
 
 Particle=Text[Uranium];
 
 R = 0.08;(*Radius of solenoid coil[m]*)
 
 l = 0.3;(*Lenght of solenoid coil[m]*)
 
 Lz = 20;(*Simulation area of longitudinal direction[m]*)
 
 clight = 2.99792458*10^8;(*Speed of light[m/s]*)
 
 m = 238*1.66053892*10^-27;(*Mass of ion[kg]*)
 
 q = 1.60217657*10^-19;(*Elementary chrge*)
 
 Inum = 34;(*Charge number*)
 
 Ene = 10*10^6;(*Kinetic energy*)
 
 gamma = q*Ene/(m*clight^2) + 1;(*Lorentz factor*)
 
 v = clight*Sqrt[1 - 1/gamma^2];(*Velosity of ion*)
 
 Xinit[1] = {0.002, 0.004};(*Initial x coodinate[m]*)
 
 Yinit[1] = {0.002, 0};(*Initial y coodinate[m]*)
 
 Pxinit[1] = {1000, 0};(*Initial vx[m/s]*)
 
 Pyinit[1] = {1000, 1000};(*Initial vy[m/s]*)
 
 BzCent[r_, z_] := NIntegrate[1/(2 \[Pi]) ((R^2 - R*r*Cos[\[Theta]])/(R^2 + r^2 - 2 R*r*Cos[\[Theta]]))*((z + l/2)/
                                Sqrt[(z + l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]] - (z - l/2)/Sqrt[(z - l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]]), {\[Theta], 0, 2 \[Pi]}];
 
 Bz0 = BzCent[0, 0];(*Scale factor*)
 
 Br[r_, z_] := NIntegrate[1/(2 \[Pi])*((R*Cos[\[Theta]])/Sqrt[(z - l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]] - (R*Cos[\[Theta]])/
                                       Sqrt[(z + l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]])/Bz0, {\[Theta], 0, 2 \[Pi]}];(*Scaled Br*)
 
 Bz[r_, z_] := NIntegrate[1/(2 \[Pi]) ((R^2 - R*r*Cos[\[Theta]])/(R^2 + r^2 - 2 R*r*Cos[\[Theta]]))*((z + l/2)/
                            Sqrt[(z + l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]] - (z - l/2)/Sqrt[(z - l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]])/Bz0, {\[Theta], 0, 2 \[Pi]}];
                            (*Scaled Bz*)
 
 Print[Text[Style["Simulation Parameters", 16, Bold, Black]]];
 (* --- information printouts *)
 headinglattice = {
     {"Ion Species",
      "Charge State[e]",
      "Kinetic Energy[eV]",
      "Lorentz Factor Gamma",
      "Lorentz Factor Beta",
      "Rigidity[Brho]",
      "Solenoid Radius[m]",
      "Solenoid Length[m]",
      "Central Solenoid Field[Tesla]",
      "Simulation Length[m]"
     },None         };
 
 outputlattice = {
     Particle,
     Inum,
     Ene,
     NumberForm[gamma,9],
     v/clight,
     m*gamma*v/q/Inum,
     R,
     l,
     1,
     Lz
     };
     
     StylePrint[ TableForm[outputlattice,TableHeadings -> headinglattice] ];
 
 Vari = {x, y, z};(*Variable of equation*)
 
 Time = {t, 0, Lz/v};(*Simulation time*)
 
 Table[Eqs[i] = {m*gamma*x''[t] == Inum*q*(y'[t]*Bz[Sqrt[x[t]^2 + y[t]^2], z[t]] - z'[t]*y[t]/Sqrt[x[t]^2 + y[t]^2]*Br[Sqrt[x[t]^2 + y[t]^2], z[t]]), 
                 m*gamma*y''[t] == Inum*q*(z'[t]*x[t]/Sqrt[x[t]^2 + y[t]^2]*Br[Sqrt[x[t]^2 + y[t]^2], z[t]] - x'[t]*Bz[Sqrt[x[t]^2 + y[t]^2], z[t]]), 
                 m*gamma*z''[t] == Inum*q*(x'[t]*y[t]/Sqrt[x[t]^2 + y[t]^2]*Br[Sqrt[x[t]^2 + y[t]^2], z[t]] - y'[t] x[t]/Sqrt[x[t]^2 + y[t]^2]*
                                   Br[Sqrt[x[t]^2 + y[t]^2], z[t]]), x[0] == Xinit[1][[i]], 
                 x'[0] == Pxinit[1][[i]], y[0] == Yinit[1][[i]], y'[0] == Pyinit[1][[i]], z[0] == -Lz/2, z'[0] == v}, {i, 1, 2}];(*Equation of motion*)
                 Table[SOL[i] = NDSolveValue[Eqs[i], Vari, Time, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 4}, MaxStepSize -> Lz/v/1000], {i, 1, 2}];
                 (*Solve equation of motion numerically*)
    
Table[xy[i] = Table[{SOL[i][[1]][t], SOL[i][[2]][t]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
Table[xz[i] = Table[{SOL[i][[3]][t], SOL[i][[1]][t]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
Table[yz[i] = Table[{SOL[i][[3]][t], SOL[i][[2]][t]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
Table[xpx[i] = Table[{SOL[i][[1]][t], SOL[i][[1]]'[t]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
Table[ypy[i] = Table[{SOL[i][[2]][t], SOL[i][[2]]'[t]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
Table[zpz[i] = Table[{SOL[i][[3]][t], SOL[i][[3]]'[t]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
    
    
Print[Text[Style["Plot of Particle Orbit in Real Space", 16, Bold, Black]]];
    
Print[XY[1] = ListPlot[Table[xy[i], {i, 1, 2}], PlotRange -> All, FrameLabel -> {"x_coordinate[m]", "y_coordinate[m]"}, 
                       AspectRatio -> 1, LabelStyle -> Directive[15], Axes -> None, PlotStyle -> PointSize[0.009], Joined -> True, Frame -> True]];(*Plot particle orbit in xy space*)
            
Print[XZ[1] = ListPlot[Table[xz[i], {i, 1, 2}], PlotRange -> All, FrameLabel -> {"z_coordinate[m]", "x_coordinate[m]"}, 
                LabelStyle -> Directive[15], PlotStyle -> PointSize[0.009], Axes -> None, Joined -> True, Frame -> True]];(*Plot particle orbit in xz space*)
            
Print[YZ[1] = ListPlot[Table[yz[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "y_coordinate[m]"}, 
                 LabelStyle -> Directive[15], PlotStyle -> PointSize[0.009], Axes -> None, Joined -> True, Frame -> True]];(*Plot particle orbit in yz space*)

Print[XPX[1] = ListPlot[Table[xpx[i], {i, 1, 2}], PlotRange -> All, FrameLabel -> {"x_coordinate[m]", "vx[m/s]"}, 
                        AspectRatio -> 1, LabelStyle -> Directive[15], PlotStyle -> PointSize[0.009], Axes -> None, Joined -> True, Frame -> True]];

Print[YPY[1] = ListPlot[Table[ypy[i], {i, 1, 2}], PlotRange -> All, FrameLabel -> {"y_coordinate[m]", "vy[m/s]"}, 
                        AspectRatio -> 1, LabelStyle -> Directive[15], PlotStyle -> PointSize[0.009], Axes -> None, Joined -> True, Frame -> True]];

Print[ZPZ[1] = ListPlot[Table[zpz[i], {i, 1, 2}], PlotRange -> All, FrameLabel -> {"z_coordinate[m]", "vz[m/s]"}, 
                        AspectRatio -> 1, LabelStyle -> Directive[15], PlotStyle -> PointSize[0.009], Axes -> None, Joined -> True, Frame -> True]];
    
Athe[r_, z_] := NIntegrate[(R^2*r)/(2 \[Pi]) ((Sin[\[Theta]]*Sin[\[Theta]])/(R^2 + r^2 - 2 R*r*Cos[\[Theta]]))*((z + l/2)/
                            Sqrt[(z + l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]] - (z - l/2)/Sqrt[(z - l/2)^2 + R^2 + r^2 - 2 R*r*Cos[\[Theta]]])/(Bz0), {\[Theta], 0, 2 \[Pi]}];
                            (*Scaled vector potential*)
    
Table[Invari[i] = Table[{SOL[i][[3]][t], m*gamma*(SOL[i][[1]][t]*SOL[i][[2]]'[t] - SOL[i][[2]][t]*SOL[i][[1]]'[t]) + 
                    Inum*q*Athe[Sqrt[SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2], SOL[i][[3]][t]]*Sqrt[SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2]}, {t,0, Lz/v, Lz/v/2000}], {i, 1, 2}];
                    (*Ptheta vs z*)

Table[MAM[i] = Table[{SOL[i][[3]][t], m*gamma*(SOL[i][[1]][t]*SOL[i][[2]]'[t] - SOL[i][[2]][t]*SOL[i][[1]]'[t])}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];
                    (*Kinetic term vs z*)
                    
Table[PT[i] = Table[{SOL[i][[3]][t], Inum*q*Athe[Sqrt[SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2], SOL[i][[3]][t]]*Sqrt[SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2]}, {t,0, Lz/v, Lz/v/2000}],{i, 1, 2}];(*Potential term vs z*)

Table[Gamma3d[i] = Table[{SOL[i][[3]][t], 1/Sqrt[1 - (SOL[i][[1]]'[t]^2 + SOL[i][[2]]'[t]^2 + SOL[i][[3]]'[t]^2)/clight^2]}, {t, 0, Lz/v, Lz/v/2000}], {i,1, 2}];(*Gamma_3dvs z*)

Table[Beta3d[i] = Table[{SOL[i][[3]][t], Sqrt[(SOL[i][[1]]'[t]^2 + SOL[i][[2]]'[t]^2 + SOL[i][[3]]'[t]^2)]/clight}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];(*Beta_3d vs z*)

Table[Gammaz[i] = Table[{SOL[i][[3]][t], 1/Sqrt[1 - (SOL[i][[3]]'[t]^2)/clight^2]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];(*Gamma_z vs z*)

Table[Betaz[i] = Table[{SOL[i][[3]][t], SOL[i][[3]]'[t]/clight}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];(*Beta_z vs z*)

Table[Rorb[i] = Table[{SOL[i][[3]][t], Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)]}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];(*r vs z*)
    
Table[Rdorb[i] = Table[{SOL[i][[3]][t], ((SOL[i][[1]][t]*SOL[i][[1]]'[t] + SOL[i][[2]][t]*SOL[i][[2]]'[t])/Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)]*1/SOL[i][[3]]'[t])}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];(*dr/dz vs z*)

Print[Text[Style["Plot of P\[Theta], Kinetic term, and Potential term", 16, Bold, Black]]];
                                         
Print[Ps[1] = ListPlot[Table[Invari[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "\!\(\*SubscriptBox[\(P\), \(\[Theta]\)]\)"}, 
                       LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];(*Plot Ptheta vs z*)

Print[Kt[1] = ListPlot[Table[MAM[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "Kitetic term"}, 
                       LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];(*Plot Knetic term vs z*)

Print[Pt[1] = ListPlot[Table[PT[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "Potential term"}, 
                       LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];(*Plot Potential term vs z*)

Print[Text[Style["Plot of Lorentz factors", 16, Bold, Black]]];                                         

Print[LG[1] = Show[ListPlot[Table[Gamma3d[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "Lorenz factor \[Gamma]"},
                            LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}],
            ListPlot[Table[Gammaz[i], {i, 1, 2}],PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "Lorenz factor \[Gamma]"}, 
                     LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]]];
                    (*Plot Lorentz gamma vs z*)
                                         
Print[LB[1] = Show[ListPlot[Table[Beta3d[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "Lorenz factor \[Beta]"}
                    , LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}],
             ListPlot[Table[Betaz[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "Lorenz factor \!\(\*SubscriptBox[\(\[Beta]\)
                    , \(b\)]\)"}, LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]]];                                         
                    (*Plot Lorentz beta vs z*)
    
Print[Text[Style["Plot of r and dr/dz orbit", 16, Bold, Black]]];
                                         
Print[Ro[1] = ListPlot[Table[Rorb[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "r_coordinate[m]"}, 
                       LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];(*Plot r vs z*)

Print[Rdo[1] = ListPlot[Table[Rdorb[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, FrameLabel -> {"z_coordinate[m]", "dr/dz"}, 
                        LabelStyle -> Directive[15], Frame -> True, Axes -> None, Joined -> True, AspectRatio -> 0.3, ImageSize -> {800, 240}]];(*Plot dr/dz vs z*)                                     
                                    
Table[RKick[i, t_] := SOL[i][[3]]'[t]*((Inum*q*Bz[Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)], SOL[i][[3]][t]])/(m*gamma)*(SOL[i][[1]][t]*SOL[i][[2]]'[t] - 
                        SOL[i][[2]][t]*SOL[i][[1]]'[t])/(Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)]*SOL[i][[3]]'[t]^2) + 
                        1/((SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)^(3/2)*SOL[i][[3]]'[t]^2)*(-(SOL[i][[1]][t]*SOL[i][[1]]'[t] + 
                        SOL[i][[2]][t]*SOL[i][[2]]'[t])^2 + (SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)*(SOL[i][[1]]'[t]^2 + 
                        SOL[i][[2]]'[t]^2)) - (Inum*q*Br[Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)], SOL[i][[3]][t]])/(m*gamma*(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)*
                        SOL[i][[3]]'[t]^3)*(SOL[i][[1]][t]*SOL[i][[1]]'[t] + SOL[i][[2]][t]*SOL[i][[2]]'[t])*(SOL[i][[2]][t]*SOL[i][[1]]'[t] - SOL[i][[1]][t]*
                        SOL[i][[2]]'[t]) - ((m*gamma*(SOL[i][[1]][0]*SOL[i][[2]]'[0] - SOL[i][[2]][0]*SOL[i][[1]]'[0]) + 
                        Inum*q*Athe[Sqrt[SOL[i][[1]][0]^2 + SOL[i][[2]][0]^2], SOL[i][[3]][0]]*Sqrt[SOL[i][[1]][0]^2 + SOL[i][[2]][0]^2])/(m*gamma*v))^2/(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)^(3/2)), {i, 1,2}];(*Radial kick in z coodrinate befor integration from simulations*)

Table[RdKick[i] =Table[{SOL[i][[3]][t],((Inum*q*Bz[Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)], SOL[i][[3]][t]])/(m*gamma)*(SOL[i][[1]][t]*SOL[i][[2]]'[t] - 
                        SOL[i][[2]][t]*SOL[i][[1]]'[t])/(Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)]*SOL[i][[3]]'[t]^2) + 
                        1/((SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)^(3/2)*SOL[i][[3]]'[t]^2)*(-(SOL[i][[1]][t]*SOL[i][[1]]'[t] + 
                        SOL[i][[2]][t]*SOL[i][[2]]'[t])^2 + (SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)*(SOL[i][[1]]'[t]^2 + 
                        SOL[i][[2]]'[t]^2)) - (Inum*q*Br[Sqrt[(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)], SOL[i][[3]][t]])/(m*gamma*(SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)*
                        SOL[i][[3]]'[t]^3)*(SOL[i][[1]][t]*SOL[i][[1]]'[t] + SOL[i][[2]][t]*SOL[i][[2]]'[t])*(SOL[i][[2]][t]*SOL[i][[1]]'[t] - SOL[i][[1]][t]*
                        SOL[i][[2]]'[t]) - ((m*gamma*(SOL[i][[1]][0]*SOL[i][[2]]'[0] - SOL[i][[2]][0]*SOL[i][[1]]'[0]) + 
                        Inum*q*Athe[Sqrt[SOL[i][[1]][0]^2 + SOL[i][[2]][0]^2], SOL[i][[3]][0]]*Sqrt[SOL[i][[1]][0]^2 + SOL[i][[2]][0]^2])/(m*gamma*v))^2/
                        (SOL[i][[1]][t]^2 + SOL[i][[2]][t]^2)^(3/2))}, {t, 0, Lz/v, Lz/v/2000}], {i, 1, 2}];    
    
Table[RDD[i] = NIntegrate[RKick[i, t], {t, 0, Lz/v}, PrecisionGoal -> 7, MaxRecursion -> 20], {i, 1, 2}];(*Radial kick in z coodrinate after integration from simulations*)
      
ClearAll[r0];

(*Analytic calculation of radial kick*)

r0[1] = Table[Sqrt[SOL[i][[1]][(Lz-3*l)/v/2]^2 + SOL[i][[2]][(Lz-3*l)/v/2]^2], {i, 1, 2}];    
    
Bz00[z_] = (z + l/2)/Sqrt[(z + l/2)^2 + R^2] - (z - l/2)/Sqrt[(z - l/2)^2 + R^2];(*Bz field*)
                                         
Bz1[z_] := ((z + l/2)/Sqrt[(z + l/2)^2 + R^2] - (z - l/2)/Sqrt[(z - l/2)^2 + R^2])/Bz00[0];(*Scaled Bz field*)
                                         
Bz2[z_] := Bz1'[z];(*Scaled Bz field*)
                                         
Bz3[z_] := Bz2'[z];(*Scaled Bz field*)
                                         
Bz4[z_] := Bz3'[z];(*Scaled Bz field*)
                                         
Bz5[z_] := Bz4'[z];(*Scaled Bz field*)
                                         
Bz6[z_] := Bz5'[z];(*Scaled Bz field*)

Brho = m*gamma*v/(Inum*q);
      
Table[RD1[i] = -1/4*NIntegrate[(Bz1[z]/Brho)^2*r0[1][[i]], {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];
                                         
Table[RD2[i] = -1/8*NIntegrate[(Bz2[z]/Brho)^2*r0[1][[i]]^3, {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];
                                         
Table[RD3[i] = -5/256*NIntegrate[(Bz3[z]/Brho)^2*r0[1][[i]]^5, {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];
                                         
Table[RD4[i] = -7/4608*NIntegrate[(Bz4[z]/Brho)^2*r0[1][[i]]^7, {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];
                                         
Table[RD5[i] = -7/98304*NIntegrate[(Bz5[z]/Brho)^2*r0[1][[i]]^9, {z, -\[Infinity], \[Infinity]}], {i, 1, 2}];
    
Print[Text[Style["Comparison of Radial kick between the theory and the simulation", 16, Bold, Black]]];                                         

Print[Rdkick[1] = ListPlot[Table[RdKick[i], {i, 1, 2}], PlotRange -> {{-5, 5}, All}, Axes -> None, Joined -> True, Frame -> True, FrameLabel -> {"z[m]","Radial Kick"}, 
                               LabelStyle -> Directive[13],  AspectRatio -> 0.75, ImageSize -> {400, 300}]]
    
    Print[Rdd[1] = Show[ListPlot[Table[{Invari[i][[1, 2]], (RD1[i] + RD2[i] + RD3[i] + RD4[i] + RD5[i])}, {i, 1, 2}], PlotRange -> {All,{-0.022,0.012}},PlotStyle -> {Blue},
                       FrameLabel -> {"\!\(\*SubscriptBox[\(P\), \(\[Theta]\)]\)", "\[CapitalDelta]r'"}, LabelStyle -> Directive[15], Frame -> True,Axes -> None],
                        ListPlot[Table[{Invari[j][[1, 2]], RDD[j]}, {j, 1, 2}], PlotRange -> All, 
                        PlotStyle -> {Red}, FrameLabel -> {"\!\(\*SubscriptBox[\(P\), \(\[Theta]\)]\)","\[CapitalDelta]r'"}, LabelStyle -> Directive[15], Frame -> True,
                                 Axes -> None]]];
                        (*Comparison of Radial kick between the theory and the simulation*)      

