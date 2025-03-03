class genlist{
	public T[] data; // Data of a generic type
	public int size => data.Length; // Length of data list
	public T this[int i] => data[i] // Creating way of indexing data
	public genlist{data = new T[0]; }
	public void add(T item){ /* add item to the list */
		T[] newdata = new T[size+1];
		System.Array.Copy(data,newdata,size);
		newdata[size]=item;
		data=newdata;
	}
}

