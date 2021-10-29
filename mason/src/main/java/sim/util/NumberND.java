/*
  Copyright 2019 by Sean Luke and George Mason University
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/

package sim.util;

/** 
    NumberND is the top-level abstract class of MASON's 2D and 3D mutable and immutable ints and doubles.
    All NumberND classes are Serializable
*/

public abstract class NumberND implements java.io.Serializable
    {
    private static final long serialVersionUID = 1;
    
    /** Returns the number of dimensions of this NumberND (normally 2 or 3) */
    public abstract int getNumDimensions();
    
    /** Returns the value at position VAL in this Number2D  (val should 0, 1, or sometimes 2) */
    public abstract double getVal(int pos);

    /** Flattens out the Number2D to an array of bytes, including the Number2D subtype. */
    public abstract byte[] toBytes();

    /** Returns whether this Number2D is mutable.  This is equivalent to asking (instanceof val MutableNumberND) */
    public boolean isMutable() { return false; }
    
    /** Returns whether this Number2D is equivalent to some other Number2D. */
    public abstract boolean equals(Object obj);

    /** Provides a hashcode for this Number2D. */
    public abstract int hashCode();
    
    /** Returns this Number2D as a nicely formatted String. */
    public abstract String toString();

    /** Returns this Number2D in mathematical coordinates as a formatted String. */
    public abstract String toCoordinates();
    
    /** Loads the value in MSB first, LSB last order. */
    protected static void intToBytes(int val, byte[] buf, int pos)
        {
        for(int i = 0; i < 4; i++)
            {
            buf[pos + (3-i)] = (byte)(val & 0xFF);
            val = val >>> 8;
            }
        }
    
    /** Loads the value in MSB first, LSB last order. */
    protected static void doubleToBytes(double val, byte[] buf, int pos)
        {
        long v = Double.doubleToRawLongBits(val);
        for(int i = 0; i < 8; i++)
            {
            buf[pos + (7-i)] = (byte)(v & 255L);
            v = v >>> 8;
            }
        }

    /** Unloads the value in MSB first, LSB last order. */
    protected static int intFromBytes(byte[] buf, int pos)
        {
        int val = 0;
        for(int i = 0; i < 4; i++)
            {
            val = val << 8;
            val = val | (buf[pos + i] & 255);
            }
        return val;
        }
    
    /** Unloads the value in MSB first, LSB last order. */
    protected static double doubleFromBytes(byte[] buf, int pos)
        {
        long val = 0;
        for(int i = 0; i < 8; i++)
            {
            val = val << 8;
            val = val | (buf[pos + i] & 255L);
            }
        return Double.longBitsToDouble(val);
        }
    
    protected static final byte TYPE_INT_2D = 0;
    protected static final byte TYPE_DOUBLE_2D = 1;
    protected static final byte TYPE_INT_3D = 2;
    protected static final byte TYPE_DOUBLE_3D = 3;
    protected static final byte TYPE_MUTABLE_INT_2D = 4;
    protected static final byte TYPE_MUTABLE_DOUBLE_2D = 5;
    protected static final byte TYPE_MUTABLE_INT_3D = 6;
    protected static final byte TYPE_MUTABLE_DOUBLE_3D = 7;

    /** Produces a Number2D from the provided array of bytes, which should have been generated by Number2D.toBytes() */
    protected static NumberND fromBytes(byte[] bytes)
        {
        switch(bytes[0])
            {
            case (TYPE_INT_2D):
                {
                return new Int2D(intFromBytes(bytes, 1), intFromBytes(bytes, 5));
                }
            case (TYPE_DOUBLE_2D):
                {
                return new Double2D(doubleFromBytes(bytes, 1), doubleFromBytes(bytes, 9));
                }
            case (TYPE_INT_3D):
                {
                return new Int3D(intFromBytes(bytes, 1), intFromBytes(bytes, 5), intFromBytes(bytes, 9));
                }
            case (TYPE_DOUBLE_3D):
                {
                return new Double3D(doubleFromBytes(bytes, 1), doubleFromBytes(bytes, 9), doubleFromBytes(bytes, 17));
                }
            case (TYPE_MUTABLE_INT_2D):
                {
                return new Int2D(intFromBytes(bytes, 1), intFromBytes(bytes, 5));
                }
            case (TYPE_MUTABLE_DOUBLE_2D):
                {
                return new Double2D(doubleFromBytes(bytes, 1), doubleFromBytes(bytes, 9));
                }
            case (TYPE_MUTABLE_INT_3D):
                {
                return new Int3D(intFromBytes(bytes, 1), intFromBytes(bytes, 5), intFromBytes(bytes, 9));
                }
            case (TYPE_MUTABLE_DOUBLE_3D):
                {
                return new Double3D(doubleFromBytes(bytes, 1), doubleFromBytes(bytes, 9), doubleFromBytes(bytes, 17));
                }
            default:
                {
                throw new RuntimeException("Invalid Number Type " + bytes[0]);
                }
            }
        }

    public static void main(String[] args)
        {
        byte[] b;

            {
            Int2D num = new Int2D(279, 321);
            NumberND num2 = fromBytes(num.toBytes());
            if (!num.equals(num2))
                System.err.println("" + num + " != " + num2);
            }

            {
            MutableInt2D num = new MutableInt2D(279, 321);
            NumberND num2 = fromBytes(num.toBytes());
            if (!num.equals(num2))
                System.err.println("" + num + " != " + num2);
            }

            {
            Int3D num = new Int3D(279, 321, -234);
            NumberND num2 = fromBytes(num.toBytes());
            if (!num.equals(num2))
                System.err.println("" + num + " != " + num2);
            }


            {
            MutableInt3D num = new MutableInt3D(279, 321, -234);
            NumberND num2 = fromBytes(num.toBytes());
            if (!num.equals(num2))
                System.err.println("" + num + " != " + num2);
            }

            {
            Double2D num = new Double2D(279.2, 321.2);
            NumberND num2 = fromBytes(num.toBytes());
            if (!num.equals(num2))
                System.err.println("" + num + " != " + num2);
            }

            {
            MutableDouble2D num = new MutableDouble2D(279.2, 321.2);
            NumberND num2 = fromBytes(num.toBytes());
            if (!num.equals(num2))
                System.err.println("" + num + " != " + num2);
            }

            {
            Double3D num = new Double3D(279.2, 321.2, -234.2);
            NumberND num2 = fromBytes(num.toBytes());
            if (!num.equals(num2))
                System.err.println("" + num + " != " + num2);
            }

            {
            MutableDouble3D num = new MutableDouble3D(279.2, 321.2, -234.2);
            NumberND num2 = fromBytes(num.toBytes());
            if (!num.equals(num2))
                System.err.println("" + num + " != " + num2);
            }

        }
    }
        