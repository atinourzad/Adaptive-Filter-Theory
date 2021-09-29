function [Ru] = RxApp(u,N)
% Ru = zeros(N,N);
% uHermitian = u';
% for i=1:N
%     for j=1:N
%         Ru(i,j) = u(i,:) * uHermitian(:,j);
%     end
% end
% Ru = Ru/N;
Ru = (u*u')/N;
end

