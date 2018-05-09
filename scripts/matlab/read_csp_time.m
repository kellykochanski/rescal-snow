function [CSP,CSP_time,H,L,D] = read_csp_time(filename)
% [CSP,CSP_time,H,L,D] = READ_CSP_TIME(FILENAME)
%    Loads a CSP-file into the arrays CSP(L,H,D) and CSP_time(L,H,D), for the state and time data respectively
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
      disp(sprintf('Model %s',char(metadata')))
    elseif (chunk == 'SIZE')
      H=fread(fid,1,'int32');
      L=fread(fid,1,'int32');
      D=fread(fid,1,'int32');
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
    CSP0=reshape(fread(fid,HLD*2,'int32'),2,L,H,D);
    CSP=reshape(CSP0(1,:,:,:),L,H,D);
    CSP_time = reshape(CSP0(2,:,:,:),L,H,D);
  end
  fclose(fid);
  return

