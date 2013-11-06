function [ A, b_L, b_U ] = dinput_encode_lin_cons( user )

Nsamples = user.Nsamples;

didx = user.idxs.d; % [ Nmodes, Nsamples + 1 ]
len_s = user.len_s;

A = zeros( Nsamples, len_s );
b_L = ones( Nsamples, 1 );
b_U = ones( Nsamples, 1 );

for i = 1:Nsamples
    A( i, didx( :, i ) ) = 1;
end