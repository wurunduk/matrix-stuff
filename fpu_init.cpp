/* ���������� ��� ����� ������ ������� */
#if defined (__i386__) || defined (__x86_64__)
void set_fpu_exception_mask (void)
{
  int cw = 0;

  /* �������� CW */
  asm volatile ("fstcw %0" : "=m" (cw) : /* no argument */);
  /* �������� ����� ���������� � SW (����� fldcw ���� ���������� �� ��� ������������ �����) */
  asm volatile ("fclex");
  /* �������� �����. ������� 6 ��� cw:
       5  4  3  2  1  0
     +--+--+--+--+--+--+
     |PM|UM|OM|ZM|DM|IM|
     +--+--+--+--+--+--+

     PM - Precision Mask
     UM - Underflow Mask
     OM - Overflow Mask
     ZM - divide-by-Zero Mask
     DM - Denormalized operand Mask
     IM - Invalid operation Mask

    Bit PM have to be leaved 1, otherwise fsin raises exception in all cases.
    Bit UM have to be leaved 1, otherwise small numbers product will
    raise exception in all cases.
    Clear lower 4 bits OM ZM DM IM.
   */
  cw &= ~0x0f;
  /* �������� CW */
  asm volatile ("fldcw %0" : /* no result */ : "m" (cw));

  /* �������� SSE capability */
  asm volatile ("mov $1, %%eax; cpuid; mov %%edx, %0" : "=m" (cw)
                : /* no argument */ : "eax", "ebx", "ecx", "edx");
  /* ��������� ��� 25 (SSE) ��� ��� 26 (SSE2) */
  if (cw & 0x6000000)
    {
      /*
       SSE - MXCSR
       Pnemonic  Bit Location    Description
         FZ      bit 15  Flush To Zero
         R+      bit 14  Round Positive
         R-      bit 13  Round Negative
         RZ      bits 13 and 14  Round To Zero
         RN      bits 13 and 14 are 0    Round To Nearest
         PM      bit 12  Precision Mask
         UM      bit 11  Underflow Mask
         OM      bit 10  Overflow Mask
         ZM      bit 9   Divide By Zero Mask
         DM      bit 8   Denormal Mask
         IM      bit 7   Invalid Operation Mask
         DAZ     bit 6   Denormals Are Zero
         PE      bit 5   Precision Flag
         UE      bit 4   Underflow Flag
         OE      bit 3   Overflow Flag
         ZE      bit 2   Divide By Zero Flag
         DE      bit 1   Denormal Flag
         IE      bit 0   Invalid Operation Flag
      */
      /* �������� MXCSR */
      asm volatile ("stmxcsr %0" : "=m" (cw) : /* no argument */);

      /* �������� ���� OM ZM DM IM. */
      cw &= ~(0x0f << 7);

      /* �������� MXCSR */
      asm volatile ("ldmxcsr %0" : /* no result */ : "m" (cw));
    }
}
#else /* ! (__i386__ || __x86_64__) */
#if defined (__powerpc__) || defined (_POWER)
void set_fpu_exception_mask (void)
{
  /* �������� ����� ���������� � FPSCR (����� ��������� ���������� ���� 
     ���������� �� ��� ������������ �����) */
  asm volatile ("mtfsb0 12");	// bit 12 (VXVC - FP invalid operation exception for invalid compare) = 0
  asm volatile ("mtfsfi 2,0");	// bits 8...11
  asm volatile ("mtfsfi 1,0");	// bits 4...7
  asm volatile ("mtfsfi 0,0");	// bits 0...3 = 0

  /* ���������� ��� ����� ������ �������
     ���� crf6:
      24 25 26 27
     +--+--+--+--+
     |VE|OE|UE|ZE|
     +--+--+--+--+

     VE - invalid operation exception enable
     OE - overflow exception enable
     UE - underflow exception enable
     ZE - zero divide exception enable

   ������������� ��� ���� ����� � �������
   asm volatile ("mtfsfi 6,15");
   �� ��� ����� ������ 3: VE OE ZE
   */
  asm volatile ("mtfsb1 24"); // bit VE
  asm volatile ("mtfsb1 25"); // bit OE
  asm volatile ("mtfsb1 27"); // bit ZE

  /* ��� XE ���������� ��������� 0, ��� ��� ����� sin
     ������ ���� ���������� */
  /* asm volatile ("mtfsb1 28"); */
}
#else /* ! __powerpc__ */
#if defined (__sparc__)
void set_fpu_exception_mask (void)
{
  int fsr = 0;

  /* �������� FSR */
  asm volatile ("st %%fsr,%0" : "=m" (fsr) : "0" (fsr));
  /* �������� ����� ���������� � FSR (����� ��������� ���������� ����
     ���������� �� ��� ������������ �����) */
  /* 0...4 bits - current exception bits,
     5...9 bits - accrued exception bits */
  fsr &= ~0x3ff;
  /* ��������� FSR */
  asm volatile ("ld %0,%%fsr" : "=m" (fsr) : "0" (fsr));

  /* �������� �����. Trap Enable Mask (FSR):
      27 26 25 24 23
     +--+--+--+--+--+
     |IM|OM|UM|ZM|XM|
     +--+--+--+--+--+

     IM - Invalid operation Mask
     OM - Overflow Mask
     UM - Underflow Mask
     ZM - divide-by-Zero Mask
     XM - ineXact result Mask

    ��� XM ���������� ��������� 1, ��� ��� ����� fsin ������ ���� ����������
    ��� UM ���������� ��������� 1, ��� ��� ����� ����� ���������� ��� ������
	� ���������� �������.
    �������� 3 ���� IM OM ZM.
   */
  fsr &= ~0xd000000;
  /* ��������� FSR */
  asm volatile ("ld %0,%%fsr" : "=m" (fsr) : "0" (fsr));
}
#else /* ! __sparc__ */
void set_fpu_exception_mask (void)
{
}
#endif /* __sparc__ */
#endif /* __powerpc__ */
#endif /* __i386__ */

