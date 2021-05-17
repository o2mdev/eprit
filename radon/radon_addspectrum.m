function [ m4d ] = radon_addspectrum( M, spectr )
%Radon_addspectrum 
%   Extends 3D phantom into 4D along a defined spectrum
%   M - matrix of an object (float, NxNxN - matrix dimension)
%   spectr - spectrum to be added
%   
%   for type Lorentz
%       - center - center of the distribution
%       - sigma - width

msize=size(M);

sigma = safeget(spectr, 'sigma', 5);
emax = msize(1);
center = safeget(spectr, 'center', emax/2);

m4d=zeros([prod(msize), msize(1)]);

x = 1:msize(1);

switch safeget(spectr, 'type', 'spectr')
    case 'Lorentzian'
        spectrum = ((1/pi)*((sigma/2)^2))./((x-center).^2+(sigma/2)^2);
    otherwise
end

A = M(:);
for ii = 1:prod(msize)
  m4d(ii,:)=spectrum * A(ii);
end

m4d = reshape(m4d, [msize, msize(1)]);
