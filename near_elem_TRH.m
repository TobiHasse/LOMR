function [ib]= near_elem_TRH(a,b)
% function by Tobias Hasse March 2015
% INPUTS:       a  vector
%               b  vector of same length as a
% OUTPUTS:      ib the index of b that is closest to the corresponding
%                  element of a
% code source:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/243878
% From: Roger Stafford
% Date: 16 Feb, 2009 02:27:01
% I interpret your question to mean that for each element of vector 'a' you
% seek the nearest element of 'b' in terms of absolute difference. The 
% following assumes that 'a' and 'b' are row vectors and that all the 
% elements of 'a' are finite. It obtains the two row vectors 'd' and 'ib' 
% of the same length as 'a'. Each element of 'd' is the absolute difference 
% between the corresponding element of 'a' and the nearest element of 'b'. 
% Each element of 'ib' is the index with respect to the 'b' vector of that 
% corresponding nearest 'b' element.    
a = a(:)';       %enforce column vectors, switch to rows
b = b(:)';
%% Row vector method
        m = size(a,2); n = size(b,2);
        [c,p] = sort([a,b]);
        q = 1:m+n; q(p) = q;
        t = cumsum(p>m);
        r = 1:n; r(t(q(m+1:m+n))) = r;
        s = t(q(1:m));
        id = r(max(s,1));
        iu = r(min(s+1,n));
        [d,it] = min([abs(a-b(id));abs(b(iu)-a)]);
        ib = id+(it-1).*(iu-id);
% %         ib = ib';

%% Column vector method  DOES NOT WORK
%         m = length(a); n = length(b);
%         [c,p] = sort([a',b']');
%         q = [1:m+n]'; q(p) = q';
%         t = cumsum(p>m);
%         r = [1:n]'; 
%         r(t(q(m+1:m+n))) = r;
%         s = t(q(1:m));
%         id = r(max(s,1));
%         iu = r(min(s+1,n));
%         [d,it] = min([abs(a-b(id));abs(b(iu)-a)]);
%         ib = id+(it-1).*(iu-id);
%         W = who;
%         putvar(W{:})
end

