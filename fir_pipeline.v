`timescale 1ns/1ps

module fir_pipeline #( parameter integer N = 100 )
( input wire clk,
  input wire rst,
  input wire vin,
  input wire signed [15:0] xin,
  output reg vout,
  output reg signed [15:0] yout
);

//______________________________________________________________________
// Helping functions
//______________________________________________________________________

function automatic signed [15:0] sat16(input signed [63:0] x);
    begin
        if (x > 64'sd32767)       sat16 = 16'sd32767;
        else if (x < -64'sd32768) sat16 = -16'sd32768;
        else                      sat16 = x[15:0];
    end
endfunction

//______________________________________________________________________

function integer clog2;
    input integer x;
    integer r;
    begin
        r = 0;
        x = x - 1;
        while (x > 0) begin
            x = x >> 1;
            r = r + 1;
        end
        clog2 = r;
    end
endfunction

//________________________________________________________________________

function integer next_pow2;
    input integer x;
    integer p;
    begin
        p = 1;
        while (p < x) p = p << 1;
        next_pow2 = p;
    end
endfunction

//____________________________________________________________________________
// Local Parameters
//____________________________________________________________________________

localparam integer P      = next_pow2(N >> 1); 
localparam integer LEVELS = clog2(P);     
localparam signed [63:0] HALF_LSB = 64'sd8192;
localparam integer       SHIFT    = 14;

//____________________________________________________________________________
// Coefficients + delay line
//____________________________________________________________________________

reg signed [15:0] xdel [0:N-1];
reg signed [15:0] h    [0:N-1];

initial begin
    $readmemh("coef_q214_hex.txt", h);
end

//_______________________________________________________________________________
// Stage 1: Update the delay line
//_______________________________________________________________________________

integer i;
reg vin1;
always @(posedge clk) begin
  if (rst) begin
    vin1 <= 1'b0;
    for (i=0; i<N; i=i+1) xdel[i] <= 16'sd0;
  end
  else begin
    vin1 <= vin;            
    if (vin) begin
      for (i=N-1; i>0; i=i-1) xdel[i] <= xdel[i-1];
      xdel[0] <= xin;
    end
  end
end

//_________________________________________________________________________________________
// Stage 2: Using symmetry property of filter coefficients 
//_________________________________________________________________________________________

reg vin2;
reg signed [16:0] xpair [0:(N >> 1) - 1];
reg signed [15:0] hmul [0:(N >> 1) - 1];

always @(posedge clk) begin
  if (rst) begin
    vin2 <= 1'b0;
    for (i = 0; i < (N >> 1); i = i+1) begin
      xpair[i] <= 17'sd0;
      hmul[i] <= 16'sd0;
    end
  end
  else begin
    vin2 <= vin1;
    for (i = 0; i < (N >> 1); i = i+1) begin
      xpair[i] <= xdel[i] + xdel[N-1-i];
      hmul[i] <= h[i]; 
    end
  end
end

//______________________________________________________________________________________________
// Stage 3: DSP input registers 
//______________________________________________________________________________________________

reg vin3;
reg signed [16:0] xmul [0:(N >> 1) - 1];
reg signed [15:0] hmul1 [0:(N >> 1) - 1];

always @(posedge clk) begin
  if (rst) begin
    vin3 <= 1'b0;
    for (i = 0; i < (N >> 1); i = i+1) begin
      xmul[i] <= 17'sd0;
      hmul1[i] <= 16'sd0;
    end
  end
  else begin
    vin3 <= vin2;
    for (i = 0; i < (N >> 1); i = i+1) begin
      xmul[i] <= xpair[i];
      hmul1[i] <= hmul[i];
    end
  end
end

//____________________________________________________________________________________________
// Stage 4: DSP input registers 
//____________________________________________________________________________________________

reg vin4;
reg signed [16:0] xmul_dsp [0:(N >> 1) - 1];
reg signed [15:0] hmul_dsp [0:(N >> 1) - 1];

always @(posedge clk) begin
  if (rst) begin
    vin4 <= 1'b0;
    for (i = 0; i < (N >> 1); i = i+1) begin
      xmul_dsp[i] <= 17'sd0;
      hmul_dsp[i] <= 16'sd0;
    end
  end
  else begin
    vin4 <= vin3;
    for (i = 0; i < (N >> 1); i = i+1) begin
      xmul_dsp[i] <= xmul[i];
      hmul_dsp[i] <= hmul1[i];
    end
  end
end

//____________________________________________________________________________________________
// Stage 5: Multiplication (DSP)
//____________________________________________________________________________________________

reg vin5;
reg signed [32:0] prod_dsp [0:(N>>1)-1];
wire signed [32:0] prod_w [0:(N>>1)-1];

genvar g;
generate
  for (g=0; g<(N>>1); g=g+1) begin : GEN_MUL
    (* use_dsp = "yes" *) 
    assign prod_w[g] = xmul_dsp[g] * hmul_dsp[g];
  end
endgenerate

always @(posedge clk) begin
  if (rst) begin
    vin5 <= 1'b0;
    for (i=0; i<(N>>1); i=i+1) prod_dsp[i] <= 33'sd0;
  end else begin
    vin5 <= vin4;
    for (i=0; i<(N>>1); i=i+1) prod_dsp[i] <= prod_w[i];
  end
end

//____________________________________________________________________________________________
// Stage 6: DSP output register 
//____________________________________________________________________________________________

reg vin6;
reg signed [32:0] prod_r [0:(N>>1)-1];

always @(posedge clk) begin
  if (rst) begin
    vin6 <= 1'b0;
    for (i=0; i<(N>>1); i=i+1) prod_r[i] <= 33'sd0;
  end else begin
    vin6 <= vin5;
    for (i=0; i<(N>>1); i=i+1) prod_r[i] <= prod_dsp[i];
  end
end

//_____________________________________________________________________________________________
// Stage 7: Sign extension and padding to next power of 2
//_____________________________________________________________________________________________

reg vin7;
reg signed [63:0] tree0 [0:P-1];

integer t;
always @(posedge clk) begin
    if (rst) begin
      vin7 <= 1'b0;
      for (t=0; t<P; t=t+1) tree0[t] <= 64'sd0;
    end
    else begin
      vin7 <= vin6;
      for (t=0; t<P; t=t+1) begin
          if (t < (N >> 1))
              tree0[t] <= {{31{prod_r[t][32]}}, prod_r[t]};
          else
              tree0[t] <= 64'sd0; 
      end
    end
end

//______________________________________________________________________________________________
// Stage 8 to 13: Adder tree
//______________________________________________________________________________________________

genvar L;
generate
    for (L=0; L<LEVELS; L=L+1) begin : GEN_LEVEL
        localparam integer OUT_SZ = (P >> (L+1));
        reg signed [63:0] treeN [0:OUT_SZ-1];
        reg vinN;

        integer k;

        if (L == 0) begin : GEN_L0
            always @(posedge clk) begin
                if (rst) begin
                  vinN <= 1'b0;
                  for (k=0; k<OUT_SZ; k=k+1)
                      treeN[k] <= 64'sd0;
                end
                else begin
                  vinN <= vin7;
                  for (k=0; k<OUT_SZ; k=k+1) begin
                      treeN[k] <= tree0[2*k] + tree0[2*k+1];
                  end
                end
            end
        end
        else begin : GEN_LGT0
            always @(posedge clk) begin
                if (rst) begin
                  vinN <= 1'b0;
                  for (k=0; k<OUT_SZ; k=k+1)
                      treeN[k] <= 64'sd0;
                end
                else begin
                  vinN <= GEN_LEVEL[L-1].vinN;
                  for (k=0; k<OUT_SZ; k=k+1) begin
                      treeN[k] <= GEN_LEVEL[L-1].treeN[2*k] + GEN_LEVEL[L-1].treeN[2*k+1];
                  end
                end
            end
        end
    end
endgenerate

//______________________________________________________________________________________________
// Stage 14:
//______________________________________________________________________________________________

reg vin_acc;
reg signed [63:0] acc_tree_r;

always @(posedge clk) begin
  if (rst) begin
    vin_acc <= 1'b0;
    acc_tree_r <= 64'sd0;
  end else begin
    vin_acc <= GEN_LEVEL[LEVELS-1].vinN;
    acc_tree_r <= GEN_LEVEL[LEVELS-1].treeN[0];
  end
end

//______________________________________________________________________________________________
// Stage 15:
//______________________________________________________________________________________________

reg vin_round1;
reg signed [63:0] acc_biased;
reg acc_sign;

always @(posedge clk) begin
  if (rst) begin
    vin_round1 <= 1'b0;
    acc_biased <= 64'sd0;
    acc_sign <= 1'b0;
  end else begin
    vin_round1 <= vin_acc;
    acc_sign <= acc_tree_r[63];
    if (acc_tree_r >= 0)
      acc_biased <= acc_tree_r + HALF_LSB;
    else
      acc_biased <= (-acc_tree_r) + HALF_LSB;
  end
end

//______________________________________________________________________________________________
// Stage 16:
//______________________________________________________________________________________________

reg vin_round2;
reg signed [63:0] y_shifted;

always @(posedge clk) begin
  if (rst) begin
    vin_round2 <= 1'b0;
    y_shifted <= 64'sd0;
  end else begin
    vin_round2 <= vin_round1;
    if (acc_sign)
      y_shifted <= -(acc_biased >>> SHIFT);
    else
      y_shifted <= acc_biased >>> SHIFT;
  end
end

//____________________________________________________________________________________________
// Stage 17:
//____________________________________________________________________________________________

always @(posedge clk) begin
  if (rst) begin
    yout <= 16'sd0;
    vout <= 1'b0;
  end else begin
    vout <= vin_round2;

    if (y_shifted > 64'sd32767)
      yout <= 16'sd32767;
    else if (y_shifted < -64'sd32768)
      yout <= -16'sd32768;
    else
      yout <= y_shifted[15:0];
  end
end

endmodule
