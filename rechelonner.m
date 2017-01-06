function y=rechelonner(x)
% ramene un signal stereo entre -1 et 1
% il faut x soit sur 2 colonnes

xmax=max(abs(x(:)));
slim = [-xmax xmax];
dx=diff(slim);

if dx==0,
    % Protect against divide-by-zero errors:
    y = zeros(size(x));
else
	y = (x-slim(1))/dx*2-1;
end
