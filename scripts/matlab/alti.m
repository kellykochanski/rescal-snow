function alti(filename)
% ALTI(FILENAME)
%    Displays the surface of a cellular space from a CSP-file.

  [csp,H,L,D]=read_csp(filename);
  
  %% sand grains
  type=0;
  b=(csp==type);
  b(:,H,:)=1;
  disp(size(b));
  
  %% number of sand grains
  nbgr=sum(sum(sum(b)));

  %% sand surface
  for k=1:L
    for i=1:D
      alt(k,i)=H-min([find(b(k,:,i))]);
    end
  end

  figure
  surf(alt)
  shading interp
  set(gca,'plotboxaspectratio',[D/H L/H 1])
  zlim([0 H])
  view(163,64)
  %view(163,80)
  camlight headlight
  %view(90,90)
  view(90,45)
  colormap(jet.^2)
  material dull
  return
