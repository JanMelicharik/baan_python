function prt_nlrm_ncp(results)
% funkce pro standardizaci vypisu vysledku posteriorni analyzy

%Vypis vysledku - beta
fprintf('Apriorni a posteriorni stredni hodnoty pro Beta (smerodatne odchylky v zavorkach)\n');
fprintf('                        Prior Beta     Posterior Beta \n');
k = length(results.bmean1);
for i=1:k
fprintf('Beta %2u             %12.4f      %12.4f        \n',[i results.b0(i) results.bmean1(i)]);
fprintf('                   (%12.4f)    (%12.4f)       \n',[results.bstd0(i) results.bstd1(i)]);
fprintf('   ------------------------------------\n');
end
fprintf('***************************************************************************\n');
fprintf('Apriorni a posteriorni stredni hodnoty pro h (smerodatne odchylky v zavorkach)\n');
fprintf('                    Prior h   Posterior h\n');
fprintf('h                    %12.4e      %12.4e \n',[results.h0  results.hmean1]);
fprintf('                    (%12.4e)    (%12.4e)\n',[results.hstd0 results.hstd1])
fprintf('***************************************************************************\n');
%    if results.pred==1
%        fprintf('                     Predikce pro x=%6.3f\n',results.xstar);
%        fprintf('E[y*]      %6.3f   \n',results.ystarmean);
%        fprintf('std[y*]    %6.3f   \n',results.ystarstd);
%    end
fprintf('Log of Marginal likelihood:    %6.3f\n\n\n',results.lmarglik);




