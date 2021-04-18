import {
  REQUEST_RECIPES,
  RECEIVE_RECIPES,
  FAILURE_RECIPES,
  ADD_RECIPE_SUCCESS,
  ADD_RECIPE_FAILURE,
  UPDATE_RECIPE_SUCCESS,
  UPDATE_RECIPE_FAILURE,
  DELETE_RECIPE_REQUEST,
  DELETE_RECIPE_SUCCESS,
  DELETE_RECIPE_FAILURE
} from "../constants/recipe";

const initialState = {
  recipeList: {
    data: [],
    isFetching: false,
    error: null
  }
};

const recipeList = (state, action) => {
  switch (action.type) {
    case REQUEST_RECIPES:
      return {
        ...state,
        isFetching: true
      };
    case RECEIVE_RECIPES:
      return {
        ...state,
        isFetching: false,
        data: action.payload
      };
    case FAILURE_RECIPES:
      return {
        ...state,
        error: action.payload
      };
    case ADD_RECIPE_SUCCESS:
      return {
        ...state,
        data: [...state.data, action.payload]
      };

    case UPDATE_RECIPE_SUCCESS:
      return {
        ...state,
        data: state.data.map(recipe =>
          recipe._id === action.payload._id ? action.payload : recipe
        )
      };

    case DELETE_RECIPE_SUCCESS:
      return {
        ...state,
        data: state.data.filter(recipe => recipe._id !== action.payload)
      };

    default:
      return state;
  }
};

export const recipeManagement = (state = initialState, action) => {
  switch (action.type) {
    default:
      return {
        recipeList: recipeList(state.recipeList, action)
      };
  }
};
