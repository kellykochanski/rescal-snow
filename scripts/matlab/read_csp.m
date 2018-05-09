function [CSP,H,L,D] = read_csp(filename)
% [CSP,H,L,D] = READ_CSP(FILENAME)
%    Loads a CSP-file into the array CSP(L,H,D).
  fprintf('loading %s\n',filename)
  fid=fopen(filename,'r');
  hd1=fread(fid,20);
  if (hd1(1:4) ~= [138;67;83;80])
    error(sprintf('%s is not in csp format',filename))
  end
  hdsize=hd1(17);
  offset=20;
  while (offset < hdsize)
    ch1=fread(fid,8);
    chksz=ch1(1);
    chunk=char(ch1(5:8)');
    mdsize=chksz-8;
    if (chunk == 'MODL')
      metadata=fread(fid,mdsize);
      fprintf('model: %s\n',char(metadata'))
    elseif (chunk == 'SIZE')
      H=fread(fid,1,'int32');
      L=fread(fid,1,'int32');
      D=fread(fid,1,'int32');
      fprintf('size: %dx%dx%d\n',H,L,D)
    elseif (chunk == 'CELL')
      szcell=fread(fid,1,'int32');
    else
      fread(fid,mdsize);
    end
    offset=offset+chksz;
  end
  HLD=H*L*D;
  csp_size=HLD*szcell;
  if (szcell == 1)
    CSP=reshape(fread(fid,csp_size),L,H,D);
  else
    CSP0=reshape(fread(fid,csp_size),szcell,L,H,D);
    CSP=reshape(CSP0(1,:,:,:),L,H,D);
  end
  fclose(fid);
  return

