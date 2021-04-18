import http from "../services/httpService";

const apiEndpoint = "http://127.0.0.1:3004/recipes";

export const recipeList = async () => {
  try {
    const { data } = await http.get(apiEndpoint);
    return data;
  } catch (error) {
    throw new Error(error);
  }
};

export const getPostById = async id => {
  try {
    const { data } = await http.get(`${apiEndpoint}/${id}`);
    return data;
  } catch (error) {
    if (error.response.status >= 400) throw new Error("Recipe not found");
  }
};

export const createRecipe = async newrecipe => {
  try {
    const { data } = await http.post(apiEndpoint, newrecipe);
    //send response to client
    return data;
  } catch (error) {
    console.log("service catch");
    throw new Error(error);
  }
};

export const updateRecipe = async (updatedrecipe, id) => {
  try {
    const { data } = await http.put(apiEndpoint + "/" + id, updatedrecipe);
    //send response to client
    return data;
  } catch (error) {
    console.log("service catch");
    throw new Error(error);
  }
};

export const deleteRecipe = async id => {
  try {
    const { data } = await http.delete(apiEndpoint + "/" + id);
    //send response to client
    return data;
  } catch (error) {
    console.log("service catch");
    throw new Error(error);
  }
};
