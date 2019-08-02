%%%coldhot
ar = [zeros(1,128) (1:78)*255/78 255*ones(1,50)]';
ag = [255*ones(1,8) (43:-1:1)*255/43 zeros(1,154) ...
      (1:43)*255/43 255*ones(1,8)]';
ab = ar(end:-1:1);

coldhot = [ar ag ab]/255;
%write_data('coldhot.bin',coldhot,'b','uchar');

%%%cold (higher is darker)
ar = zeros(256,1);
ag = [255*ones(1,16) (86:-1:1)*255/86 zeros(1,154)]';
ab = flipud([zeros(1,128) (1:78)*255/78 255*ones(1,50)]');

cold = [ar ag ab]/255;
%write_data('cold.bin',cold,'b','uchar');

%%%cold (higher is lighter)
ar = zeros(256,1);
ag = [zeros(1,154) (1:86)*255/86 255*ones(1,16)]';
ab = 50 + [zeros(1,128) (1:78)*204/78 204*ones(1,50)]';ab(1) = 0;

cold = [ar ag ab]/255;
%write_data('cold2.bin',cold,'b','uchar');

%%%hot (note: also defined by matlab...with white at top instead of yellow)
ar = [(0:153)*255/154 ones(1,102)*255]'/255;
ag = [zeros(1,154) (1:86)*255/86 255*ones(1,16)]'/255;
ab = [zeros(1,256)]';

hot = [ar ag ab];
%write_data('hot.bin',coldhot,'b','uchar');