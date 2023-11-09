function[E,psi]=Sch1D(z,V0,Mass,n,kvec,super_lattice_N)

  h=6.62606896E-34;
  hbar=h/(2*pi);
  e=1.602176487E-19;
  m0=9.10938188E-31;

  Nz=length(z);
  dz=z(2)-z(1);

  DZ2 =(-2)*diag(ones(1,Nz)) + (1)*diag(ones(1,Nz-1),-1) + (1)*diag(ones(1,Nz-1),1);
  period = Nz*dz*(super_lattice_N+1);
  DZ2(1,Nz) = exp(-1i*kvec*period); DZ2(Nz,1) = exp(1i*kvec*period);
  DZ2=DZ2/dz^2;

  H = -hbar^2/(2*m0*Mass) * DZ2  +  diag(V0*e);

  H = sparse(H);
  [psi,Energy] = eigs(H,n,'SM');
  E = diag(Energy)/e ;
  E=real(E);

  for i=1:n
      psi(:,i)=psi(:,i)/sqrt(trapz(z',abs(psi(:,i)).^2));
  end

  if E(1)>E(2)
    psi=psi(:,end:-1:1);
    E=E(end:-1:1);
  end

end
