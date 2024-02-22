;+
; :Description:
;    This code is used to compute the distance codingdiversity (dcd).
;
; :Params:
;    s_win: s_win is the window size
;    bin: bin is the number of gray levels
;
; :Keywords:
;    Simpson: the default is to compute the Shannon-Wiener index
;    
; dcd calculating
function dcd_cal, data, length_win, Simpson=Simpson
  h=histogram(data, binsize=1, reverse_indices=r, /nan)
  index=where(h gt 0)  
  pi=!null  

  foreach i, index do begin
    
    sub=r[r[i]:r[i+1]-1]  
    code = shift(sub, -1) - sub 
    code[-1]=code[-1]+length_win 

    h1=histogram(code, binsize=1)
    index1=where(h1 gt 0)
    pi=[pi,h1[index1]]
  endforeach

  pi=pi/float(length_win)
  if keyword_set(Simpson) then return, 1.0-total(pi*pi) $ 
  else return, total(-pi*alog(pi)) 

end

; circle arrangement
function circle, data, s_win
  A=reform(data,[1,s_win*s_win])
  CASE (s_win) OF
    3: A=[A[4], A[5], A[8], A[7], A[6], A[3], a[0], A[1],A[2]]
    5: A=[A[12], A[13], A[14], A[19], A[18], A[24], A[23], A[17],A[22],A[21], A[16], A[20], A[15], A[11], A[10], A[5], A[6], A[0], A[1], A[7], A[2], A[3], A[8], A[4],A[9]]
    7: A=[A[24], A[25], A[26], A[27], A[34], A[33], A[41], A[48],A[40],A[32], A[47], A[39], A[46], A[45], A[38], A[31], A[44],A[37],A[43], A[42], A[36], A[30], A[35], A[29],A[28], A[21], A[22], A[23], A[14], A[7], A[15], A[16],A[8],A[0], A[1], A[9], A[2], A[3], A[10], A[17], A[4],A[11],A[5], A[6], A[12], A[18], A[13], A[19], A[20]]
    ELSE: BEGIN
    END
  ENDCASE
  return, A
end

;serpentine arrangement
function serpentine, data
  sz=size(data)
  n=sz(1)
  A=fltarr(n,n)
  for i=0,n-1 do begin
    if i mod 2 eq 0 then begin
      A[*,i]=data[*,i]
    endif else begin
      A[*,i]=reverse(data[*,i])
    endelse
  endfor
  A=reform(A,[1,n*n])
  return, A
end

; Call the function to compute the dcd of the whole image
pro dcd
  s_win=5
  bin=12
  length_win=s_win*s_win
  fn_read=dialog_pickfile(title='data')
  img=read_image(fn_read)
  img = read_tiff('E:\application_test\cropland\test\test_crop.tif', GEOTIFF=GeoKeys)
  img=mean(img,dimension = 1)
  a=max(img)
  b=min(img)
  img=ceil(img/((a-b)/bin))
  sz=size(img)
  n1=sz[1]
  n2=sz[2]
  help,n1,n2
  dcd_img=fltarr(n1,n2)
  half=(s_win-1)/2
  help,half
  for i=half,n1-half-1 do begin

    for j=half,n2-half-1 do begin
      S=dcd_cal(img[i-half:i+half,j-half:j+half], length_win)
      C=dcd_cal(transpose(img[i-half:i+half,j-half:j+half]), length_win)
      SS=dcd_cal(serpentine(img[i-half:i+half,j-half:j+half]),length_win)
      CS=dcd_cal(serpentine(transpose(img[i-half:i+half,j-half:j+half])), length_win)
      CIRCLE_i=dcd_cal(circle(img[i-half:i+half,j-half:j+half],s_win),length_win)
      dcd_img[i,j]=(S+C+SS+CS+CIRCLE_i)/float(5)
      print,dcd_img[i,j]
    endfor
    
  endfor
  
  fn_write=dialog_pickfile(title='dcd_img')
  write_image, fn_write,'TIFF',dcd_img
 
end

