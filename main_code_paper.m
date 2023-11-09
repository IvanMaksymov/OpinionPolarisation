######################################
# This Matlab/Octave code produces one
# of the figures in the manuscript
# entitled "A Quantum-Mechanical Approach
# to Asymmetric Opinion Polarisation".
# The other figures can be easily obtained
# by substituting the model parameters
# presented in the main text of the article.
# The code calls the function Sch1D
# that implements the finite-difference
# model presented in the article. Note
# that 20 wave vector points are used
# in this code compared with 50 points
# used in the article, which is done to
# reduce the computational effort.
# Using a workstation computer, you
# should expect a 10-15 minutes
# run if you want to use 50 points.
# Octave uses # to comment the code.
# Matlab users should replace it for %.
######################################

clc; close all; clear all;

fnt_size = 22;

dz=2e-10;#discretisation step for the finite-difference algorithm, in metres
n=30;
mass = 0.067;
barrier=0.1;#the depth of the potential wells in eV

super_lattice_N = 10;#the number of wells in the superlattice
periodic_structure = [20 20 20 20 20 20 20 20 20 20 19.999 20 20 20 20 20 20 20 20 20 20]*1e-9;#the sequence of potential well
wells = [0.7 0 0.7 0 0.7 0 0.7 0 0.7 0 0.6999 0 0.7 0 0.7 0 0.7 0 0.7 0 0.7]*barrier;#the depth of the potential wells

for i=1:length(periodic_structure)
  if i==1
    zv{1} = 0:dz:periodic_structure(1);
    z=zv{1};
    V0=zv{1}*0+wells(1);
  else
    zv{i} = (z(end)+dz):dz:(z(end) + periodic_structure(i));
    z  = [ z  zv{i} ];
    V0 = [ V0   (zv{i}*0+1) * wells(i)  ];
  end
end

f = figure(1);

Nz=length(z);  dz=z(2)-z(1);
Nkv = 20; period = Nz*dz*(super_lattice_N+1);#Nkv is the number of wavevector points
for ii=-Nkv+1:Nkv

  kvec = ((ii-1)/Nkv)*(pi/period);

  [E1,psi1] = Sch1D(z,V0,mass,n,kvec,super_lattice_N);

  E=nan(n,1);
  E(1:length(E1),1)=E1;

  subplot(1,9,[1]); hold on;

  for jj=1:length(E1)
    line([0 1],[E1(jj) E1(jj)],'linewidth',1,'color','r');
  end;

  ylim([min(V0)-0.0001 max(V0)+0.1]);
  xlim([0 1]);
  set(gca,'linewidth',2);
  set(gca,'fontsize',fnt_size);
  set(gca,'xtick',[]);

  a = axis;
  plot([a(1) a(2)],[a(3) a(3)],'w','linewidth',2);

  subplot(1,9,[2 4 5 6]); hold on;

  for jj=1:length(E1)
    line([0+(ii-1) 1+(ii-1)],[E1(jj) E1(jj)],'linewidth',2,'color','r');
  end;
  ylim([min(V0)-0.0001 max(V0)+0.1]);

  xlim([-21 21]);
  xlabel('$kaN_{cell}/\pi$', "interpreter", "latex");
  set(gca,'linewidth',2);
  set(gca,'fontsize',fnt_size);
  box on;
  xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50]);
  xticklabels({'-1.0','-0.8','-0.6','-0.4','-0.2','0.0','0.2', '0.4', '0.6', '0.8', '1.0'});

end;


##################### DEFECT (opposing view) ####################

super_lattice_N = 10;
periodic_structure = [20 20 20 20 20 20 20 20 20 20 9.999 20 20 20 20 20 20 20 20 20 20]*1e-9;
wells = [0.7 0 0.7 0 0.7 0 0.7 0 0.7 0 1 0 0.7 0 0.7 0 0.7 0 0.7 0 0.7]*barrier;

for i=1:length(periodic_structure)
  if i==1
    zv{1} = 0:dz:periodic_structure(1);
    z=zv{1};
    V0=zv{1}*0+wells(1);
  else
    zv{i} = (z(end)+dz):dz:(z(end) + periodic_structure(i));
    z  = [ z  zv{i} ];
    V0 = [ V0   (zv{i}*0+1) * wells(i)  ];
  end
end

Nz=length(z);  dz=z(2)-z(1);
Nkv = 20; period = Nz*dz*(super_lattice_N+1);
for ii=-Nkv+1:Nkv

  kvec = ((ii-1)/Nkv)*(pi/period);

  [E1,psi1] = Sch1D(z,V0,mass,n,kvec,super_lattice_N);

  E=nan(n,1);
  E(1:length(E1),1)=E1;

  subplot(1,9,[1]); hold on;

  for jj=1:length(E1)
    line([1 2],[E1(jj) E1(jj)],'linewidth',1,'color','b','linestyle','-');
  end;

  ylim([-0.0001 0.1]);
  xlim([0 2]);

  set(gca,'linewidth',2);
  set(gca,'fontsize',fnt_size);
  set(gca,'xtick',[]);

  a = axis;
  plot([a(1) a(2)],[a(3) a(3)],'w','linewidth',2);

  subplot(1,9,[2 4 5 6]); hold on;

  for jj=1:length(E1)
    line([0+(ii-1) 1+(ii-1)],[E1(jj) E1(jj)],'linewidth',2,'color','b','linestyle',':');
  end;
  ylim([-0.0001 0.1]);

  xlim([-21 21]);
  xlabel('$kaN_{cell}/\pi$', "interpreter", "latex");
  set(gca,'linewidth',2);
  set(gca,'fontsize',fnt_size);
  box on;
  xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50]);
  xticklabels({'-1.0','-0.8','-0.6','-0.4','-0.2','0.0','0.2', '0.4', '0.6', '0.8', '1.0'});

end;

figure(f,"position",get(0,"screensize"))
set(0,'DefaultAxesTitleFontWeight','normal');

set(f,'PaperSize',[20 15]);
print (f, "Figure_test.pdf", "-dpdf");

stop


